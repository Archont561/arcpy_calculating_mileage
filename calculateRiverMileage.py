import arcpy, math, json, random
from collections import defaultdict
from pathlib import Path

WORKSPACE = Path(__file__).parent / "temp"


def calculate_distance(x1: float, y1: float, x2: float, y2: float) -> float:
    y_delta = y2 - y1
    x_delta = x2 - x1
    distance = math.sqrt(y_delta**2 + x_delta**2)
    return distance


def order_lines(line_features: str | arcpy.Result) -> None:
    line_features = f"{str(line_features).split('.')[0]}"
    lines_dangles = arcpy.management.FeatureVerticesToPoints(
        in_features=line_features,
        out_feature_class=f"{line_features}_line_dangles",
        point_location="DANGLE"
    )
    with arcpy.da.SearchCursor(lines_dangles, ["OID@", "ORIG_FID"]) as cursor:
        temp = defaultdict(list)
        for fid, line_fid in cursor: temp[line_fid].append(str(fid))
        double_dangle_lines_ids = []
        for line_fid, fids in temp.items(): 
            if len(fids) == 2: double_dangle_lines_ids.append(str(line_fid))
    
    lines_beginings = arcpy.management.FeatureVerticesToPoints(
        in_features=line_features,
        out_feature_class=f"{line_features}_line_beginings",
        point_location="START"
    )
    lines_dangling_begingings = arcpy.management.SelectLayerByLocation(
        in_layer=lines_dangles,
        overlap_type="ARE_IDENTICAL_TO",
        select_features=lines_beginings,
        search_distance="0 Meters",
        selection_type="NEW_SELECTION",
    )
    with arcpy.da.SearchCursor(lines_dangling_begingings, ["ORIG_FID"], 
        where_clause=f"ORIG_FID NOT IN ({','.join(double_dangle_lines_ids)})") as cursor:
        lines_to_reverse_ids = [str(row[0]) for row in cursor]

    if lines_to_reverse_ids:
        lines_to_reverse = arcpy.management.SelectLayerByAttribute(
            in_layer_or_view=line_features,
            selection_type="NEW_SELECTION",
            where_clause=f"FID IN ({','.join(lines_to_reverse_ids)})",
        )
        arcpy.edit.FlipLine(lines_to_reverse)

    for path in [lines_dangles, lines_beginings]: arcpy.management.Delete(path)

    arcpy.management.SelectLayerByAttribute(
        in_layer_or_view=line_features,
        selection_type="CLEAR_SELECTION",
        where_clause="1=1",
    )


def rank_rivers(rivernet: str | arcpy.Result, rank_field_name: str) -> None:
    river_dangles = arcpy.management.FeatureVerticesToPoints(
        in_features=rivernet,
        out_feature_class=f"{rivernet.split('.')[0]}_dangling_vertices",
        point_location="DANGLE"
    )
    with arcpy.da.SearchCursor(river_dangles, ["OID@", "ORIG_FID"]) as cursor:
        temp = defaultdict(list)
        for fid, river_fid in cursor: temp[river_fid].append(str(fid))
        main_rivers_ids = []
        for river_fid, fids in temp.items(): 
            if len(fids) == 2: main_rivers_ids.append(str(river_fid))

    arcpy.management.Delete(river_dangles)

    main_rivers_selection = arcpy.management.SelectLayerByAttribute(
        in_layer_or_view=rivernet,
        selection_type="NEW_SELECTION",
        where_clause=f"FID IN ({','.join(main_rivers_ids)})",
    )

    rank = 1
    arcpy.management.CalculateField(
        in_table=main_rivers_selection,
        field=rank_field_name,
        field_type="SHORT",
        expression=rank,
        expression_type="PYTHON3"
    )

    tributaries_count = 1
    tributaries_selection = None
    while tributaries_count > 0:
        tributaries_selection = arcpy.management.SelectLayerByLocation(
            in_layer=rivernet,
            overlap_type="INTERSECT",
            select_features=tributaries_selection or main_rivers_selection,
            search_distance="0 Meters",
            selection_type="NEW_SELECTION"
        )
        tributaries_selection = arcpy.management.SelectLayerByAttribute(
            in_layer_or_view=tributaries_selection,
            selection_type="SUBSET_SELECTION",
            where_clause=f"{rank_field_name} = 0"
        )
        tributaries_count = int(arcpy.management.GetCount(tributaries_selection)[0])
        rank += 1
        arcpy.management.CalculateField(
            in_table=tributaries_selection,
            field=rank_field_name,
            expression=rank,
            expression_type="PYTHON3"
        )

    arcpy.management.SelectLayerByAttribute(
        in_layer_or_view=rivernet,
        selection_type="CLEAR_SELECTION",
        where_clause="1=1"
    )

 
def get_intersection_points(line_feature: str | arcpy.Result) -> arcpy.Result:
    start_end_points = arcpy.management.FeatureVerticesToPoints(
        in_features=line_feature,
        out_feature_class="start_end_points",
        point_location="BOTH_ENDS"
    )
    dangling_points = arcpy.management.FeatureVerticesToPoints(
        in_features=line_feature,
        out_feature_class="dangling_points",
        point_location="DANGLE"
    )
    intersection_points = arcpy.analysis.SymDiff(
        in_features=start_end_points,
        update_features=dangling_points,
        out_feature_class=f"{str(line_feature).split('.')[0]}_intersection",
        join_attributes="ONLY_FID"
    )
    arcpy.management.Delete(start_end_points)
    arcpy.management.Delete(dangling_points)
    return intersection_points


def distance_from_begining(
    lines: str | arcpy.Result, 
    points: str | arcpy.Result,
    distance_field: str = "Distance",
    tolerance: float = 0.001,
) -> arcpy.Result:
    extension = Path(arcpy.env.workspace).suffix
    distances_table = arcpy.management.CreateTable(
        out_name=f"distance_table" + ("" if extension.endswith("gdb") else ".dbf"))
    arcpy.management.AddFields(
        in_table=distances_table,
        field_description=[
            (distance_field, "DOUBLE"),
            ((line_id_field := "LINE_FID"), "SHORT"),
            ((point_id_field := "POINT_FID"), "SHORT")
        ]
    )
    arcpy.management.DeleteField(
        in_table=distances_table, 
        drop_field=[distance_field, line_id_field, point_id_field], 
        method="KEEP_FIELDS"
    )
    with (arcpy.da.SearchCursor(lines, ["OID@", "SHAPE@"]) as line_cursor,
        arcpy.da.InsertCursor(distances_table, [line_id_field, point_id_field, distance_field]) as table_cursor):
        for line_id, line_shape in line_cursor:
            current_line = arcpy.management.SelectLayerByAttribute(
                in_layer_or_view=lines,
                selection_type="NEW_SELECTION",
                where_clause=f"FID = {line_id}",
            )
            current_intersections = arcpy.management.SelectLayerByLocation(
                in_layer=points,
                select_features=current_line,
                overlap_type="INTERSECT",
                selection_type="NEW_SELECTION",
                search_distance=f"{tolerance} Meters",
            )
            with arcpy.da.SearchCursor(current_intersections, ["OID@", "SHAPE@"]) as point_cursor:
                for point_id, point_shape in point_cursor:
                    table_cursor.insertRow([line_id, point_id, line_shape.measureOnLine(point_shape)])
    return distances_table


def calculate_rivernet_mileage(
    rivernet: str | arcpy.Result,
    features_to_locate: str | arcpy.Result, 
    pour_points: str | arcpy.Result,
    output_name: str = "",
    rank_field: str = "",
    tolerance: int = 0.001,
) -> arcpy.Result:

    if not rank_field:
        rank_field = "Rank"
        rank_rivers(rivernet, rank_field)

    with arcpy.da.SearchCursor(rivernet, ["OID@", rank_field]) as river_cursor:
        river_fid_by_ranks = defaultdict(list)
        for fid, rank in river_cursor: river_fid_by_ranks[rank].append(str(fid))
        minimum_rank = min(river_fid_by_ranks.keys())
        main_rivers_ids = []
        for rank, fids in river_fid_by_ranks.items():
             if rank == minimum_rank: main_rivers_ids.extend(fids)

    main_rivers = arcpy.management.SelectLayerByAttribute(
        in_layer_or_view=rivernet,
        selection_type="NEW_SELECTION",
        where_clause=f"FID IN ({','.join(main_rivers_ids)})"
    )
    with arcpy.da.SearchCursor(main_rivers, ["OID@", "SHAPE@"]) as main_river_cursor:
        main_rivers_to_reverse_ids = []
        for fid, shape in main_river_cursor:
            current_main_river = arcpy.management.SelectLayerByAttribute(
                in_layer_or_view=rivernet,
                selection_type="NEW_SELECTION",
                where_clause=f"FID = {fid}",
            )
            current_pour_point = arcpy.management.SelectLayerByLocation(
                in_layer=pour_points,
                select_features=current_main_river,
                overlap_type="INTERSECT",
                selection_type="NEW_SELECTION",
                search_distance="0 Meters",
            )
            with arcpy.da.SearchCursor(current_pour_point, ["SHAPE@XY"]) as point_cursor:
                try: point_shape = next(point_cursor)[0]
                except StopIteration as e: continue
            
            river_first_point = (shape.firstPoint.X, shape.firstPoint.Y)
            river_last_point = (shape.lastPoint.X, shape.lastPoint.Y)
            pour_point_river_start_point_distance = calculate_distance(*point_shape, *river_first_point)
            pour_point_river_last_point_distance = calculate_distance(*point_shape, *river_last_point) 
            if pour_point_river_start_point_distance > pour_point_river_last_point_distance:
                main_rivers_to_reverse_ids.append(str(fid))
                
    if main_rivers_to_reverse_ids:
        main_rivers_to_reverse = arcpy.management.SelectLayerByAttribute(
            in_layer_or_view=rivernet,
            selection_type="NEW_SELECTION",
            where_clause=f"FID IN ({','.join(main_rivers_to_reverse_ids)})"
        )
        arcpy.edit.FlipLine(main_rivers_to_reverse)
    
    order_lines(rivernet)
    distance_field = "Distance"
    intersections = get_intersection_points(rivernet)
    intersections_mileages = distance_from_begining(rivernet, intersections, distance_field)
    with arcpy.da.UpdateCursor(intersections_mileages, 
        ["LINE_FID", "POINT_FID", distance_field]) as intersections_cursor:
        inter_mileages = dict()
        for river_fid, point_fid, distance in intersections_cursor:
            if point_fid not in inter_mileages:
                inter_mileages[point_fid] = {"main": river_fid, "distance": distance}
            else:
                if distance > inter_mileages[point_fid]["distance"]:
                    inter_mileages[point_fid] = {
                        "main": river_fid, 
                        "tributary": inter_mileages[point_fid]["main"], 
                        "distance": distance
                    }
                else: inter_mileages[point_fid]["tributary"] = river_fid
        inter_mileages = [(v["distance"], v["main"], v["tributary"]) for v in inter_mileages.values()]

    for path in [intersections, intersections_mileages]: 
        arcpy.management.Delete(path)

    features_to_locate_mileages = distance_from_begining(
        rivernet, features_to_locate, distance_field, tolerance=tolerance)    

    full_distance_field = "FullDist"
    arcpy.management.CalculateField(
        in_table=features_to_locate_mileages,
        field=full_distance_field,
        field_type="DOUBLE",
        expression=0,
        expression_type="PYTHON3"
    )
    path_field = "Path"
    arcpy.management.CalculateField(
        in_table=features_to_locate_mileages,
        field=path_field,
        field_type="TEXT",
        expression="''",
        expression_type="PYTHON3"
    )

    with arcpy.da.UpdateCursor(features_to_locate_mileages,
        ["OID@", distance_field, "LINE_FID", full_distance_field, path_field],
        where_clause=f"{distance_field} > 0") as points_cursor:
        for fid, distance, river_id, _, _ in points_cursor:
            next_river_rank = float("inf")
            full_distance = distance
            current_river_id = river_id
            current_path = [str(river_id)]
            while True:
                current_row = [row for row in inter_mileages if row[-1] == current_river_id]
                if not current_row: break
                inter_distance, next_river_id, current_river_id = current_row[0]
                full_distance += inter_distance
                current_path.append(str(next_river_id))
                current_river_id = next_river_id

            points_cursor.updateRow([fid, distance, river_id, full_distance, ",".join(current_path)])

    return features_to_locate_mileages


def points_from_mileage(
    mileage: str | arcpy.Result, 
    rivernet: str | arcpy.Result, 
    river_id_field: str,
    distance_field: str
) -> arcpy.Result:
    mileage_points = arcpy.management.CreateFeatureclass(
        out_name="mileage_points",
        geometry_type="POINT",
        spatial_reference=arcpy.Describe(rivernet).spatialReference,
    )
    arcpy.management.AddField(
        in_table=mileage_points,
        field_name=river_id_field,
        field_type="SHORT",
    )
    arcpy.management.AddField(
        in_table=mileage_points,
        field_name=distance_field,
        field_type="DOUBLE",
    )
    with (arcpy.da.SearchCursor(rivernet, ["OID@", "SHAPE@"]) as current_river_cursor,
        arcpy.da.SearchCursor(mileage, [river_id_field, distance_field]) as mileage_cursor,
        arcpy.da.InsertCursor(mileage_points, [river_id_field, "SHAPE@", distance_field]) as mileage_points_cursor):
        river_distances = [(river_id, distance) for river_id, distance in mileage_cursor]
        for river_id, river_shape in current_river_cursor:
            current_river_distances = [distance for rid, distance in river_distances if river_id == rid]
            for distance in current_river_distances:
                mileage_points_cursor.insertRow(
                    [river_id, river_shape.positionAlongLine(distance), distance])
    
    arcpy.management.DeleteField(
        in_table=mileage_points,
        drop_field=['Id']
    )
    return mileage_points


def generate_points_along_line(
    line_features: str | arcpy.Result, 
    points_per_feature: int
) -> arcpy.Result:
    random_line_points = arcpy.management.CreateFeatureclass(
        out_name="random_points",
        geometry_type="POINT",
        spatial_reference=arcpy.Describe(line_features).spatialReference,
    )
    with (arcpy.da.InsertCursor(random_line_points, ["SHAPE@"]) as random_points_cursor,
        arcpy.da.SearchCursor(line_features, ["SHAPE@"]) as line_features_cursor):
        for [line_feature] in line_features_cursor:
            for _ in range(points_per_feature):
                random_points_cursor.insertRow(
                    line_feature.positionAlongLine(random.random()*line_feature.length))
    return random_line_points


def main():
    with arcpy.EnvManager(workspace=str(WORKSPACE), overwriteOutput=True):
        # output_name = calculate_rivernet_mileage(
        #     rivernet="rivers.shp",
        #     features_to_locate="random_points.shp",
        #     pour_points="pour_points.shp",
        #     rank_field="Rank",
        # )
        points_from_mileage("distance_table.dbf", "rivers.shp", "LINE_FID", "Distance")
        print("END")


if __name__ == "__main__":
    main()