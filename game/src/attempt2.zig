// package main

// import (
//     "fmt"
//     "sort"
// )

// type (
//     Point [3]float64
//     Face  []int

//     Edge struct {
//         point1 int   // point number 1
//         point2 int   // point number 2
//         face1 int   // face number 1
//         face2 int   // face number 2
//         centerPoint  Point // center point
//     }

//     PointEx struct {
//         p Point
//         n int
//     }
// )

// func sumPoint(p1, p2 Point) Point {
//     sp := Point{}
//     for i := 0; i < 3; i++ {
//         sp[i] = p1[i] + p2[i]
//     }
//     return sp
// }

// func mulPoint(p Point, m float64) Point {
//     mp := Point{}
//     for i := 0; i < 3; i++ {
//         mp[i] = p[i] * m
//     }
//     return mp
// }

// func divPoint(p Point, d float64) Point {
//     return mulPoint(p, 1.0/d)
// }

// func centerPoint(p1, p2 Point) Point {
//     return divPoint(sumPoint(p1, p2), 2)
// }

// func getFacePoints(inputPoints []Point, inputFaces []Face) []Point {
//     facePoints := make([]Point, len(inputFaces))
//     for i, currFace := range inputFaces {
//         facePoint := Point{}
//         for _, centerPointi := range currFace {
//             currPoint := inputPoints[centerPointi]
//             facePoint = sumPoint(facePoint, currPoint)
//         }
//         facePoint = divPoint(facePoint, float64(len(currFace)))
//         facePoints[i] = facePoint
//     }
//     return facePoints
// }

// func getEdgesFaces(inputPoints []Point, inputFaces []Face) []Edge {
//     var edges [][3]int
//     for faceNum, face := range inputFaces {
//         numPoints := len(face)
//         for pointIndex := 0; pointIndex < numPoints; pointIndex++ {
//             point1 := face[pointIndex]
//             var point2 int
//             if pointIndex < numPoints-1 {
//                 point2 = face[pointIndex+1]
//             } else {
//                 point2 = face[0]
//             }
//             if point1 > point2 {
//                 point1, point2 = point2, point1
//             }
//             edges = append(edges, [3]int{point1, point2, faceNum})
//         }
//     }
//     sort.Slice(edges, func(i, j int) bool {
//         if edges[i][0] == edges[j][0] {
//             if edges[i][1] == edges[j][1] {
//                 return edges[i][2] < edges[j][2]
//             }
//             return edges[i][1] < edges[j][1]
//         }
//         return edges[i][0] < edges[j][0]
//     })
//     numEdges := len(edges)
//     eIndex := 0
//     var mergedEdges [][4]int
//     for eIndex < numEdges {
//         e1 := edges[eIndex]
//         if eIndex < numEdges-1 {
//             e2 := edges[eIndex+1]
//             if e1[0] == e2[0] && e1[1] == e2[1] {
//                 mergedEdges = append(mergedEdges, [4]int{e1[0], e1[1], e1[2], e2[2]})
//                 eIndex += 2
//             } else {
//                 mergedEdges = append(mergedEdges, [4]int{e1[0], e1[1], e1[2], -1})
//                 eIndex++
//             }
//         } else {
//             mergedEdges = append(mergedEdges, [4]int{e1[0], e1[1], e1[2], -1})
//             eIndex++
//         }
//     }
//     var edgesCenters []Edge
//     for _, me := range mergedEdges {
//         p1 := inputPoints[me[0]]
//         p2 := inputPoints[me[1]]
//         centerPoint := centerPoint(p1, p2)
//         edgesCenters = append(edgesCenters, Edge{me[0], me[1], me[2], me[3], centerPoint})
//     }
//     return edgesCenters
// }

// func getEdgePoints(inputPoints []Point, edgesFaces []Edge, facePoints []Point) []Point {
//     edgePoints := make([]Point, len(edgesFaces))
//     for i, edge := range edgesFaces {
//         centerPoint := edge.centerPoint
//         fp1 := facePoints[edge.face1]
//         var fp2 Point
//         if edge.face2 == -1 {
//             fp2 = fp1
//         } else {
//             fp2 = facePoints[edge.face2]
//         }
//         cfp := centerPoint(fp1, fp2)
//         edgePoints[i] = centerPoint(centerPoint, cfp)
//     }
//     return edgePoints
// }

// func getAvgFacePoints(inputPoints []Point, inputFaces []Face, facePoints []Point) []Point {
//     numPoints := len(inputPoints)
//     tempPoints := make([]PointEx, numPoints)
//     for faceNum := range inputFaces {
//         fp := facePoints[faceNum]
//         for _, pointNum := range inputFaces[faceNum] {
//             tp := tempPoints[pointNum].p
//             tempPoints[pointNum].p = sumPoint(tp, fp)
//             tempPoints[pointNum].n++
//         }
//     }
//     avgFacePoints := make([]Point, numPoints)
//     for i, tp := range tempPoints {
//         avgFacePoints[i] = divPoint(tp.p, float64(tp.n))
//     }
//     return avgFacePoints
// }

// func getAvgMidEdges(inputPoints []Point, edgesFaces []Edge) []Point {
//     numPoints := len(inputPoints)
//     tempPoints := make([]PointEx, numPoints)
//     for _, edge := range edgesFaces {
//         centerPoint := edge.centerPoint
//         for _, pointNum := range []int{edge.point1, edge.point2} {
//             tp := tempPoints[pointNum].p
//             tempPoints[pointNum].p = sumPoint(tp, centerPoint)
//             tempPoints[pointNum].n++
//         }
//     }
//     avgMidEdges := make([]Point, len(tempPoints))
//     for i, tp := range tempPoints {
//         avgMidEdges[i] = divPoint(tp.p, float64(tp.n))
//     }
//     return avgMidEdges
// }

// func getPointsFaces(inputPoints []Point, inputFaces []Face) []int {
//     numPoints := len(inputPoints)
//     pointsFaces := make([]int, numPoints)
//     for faceNum := range inputFaces {
//         for _, pointNum := range inputFaces[faceNum] {
//             pointsFaces[pointNum]++
//         }
//     }
//     return pointsFaces
// }

// func getNewPoints(inputPoints []Point, pointsFaces []int, avgFacePoints, avgMidEdges []Point) []Point {
//     newPoints := make([]Point, len(inputPoints))
//     for pointNum := range inputPoints {
//         n := float64(pointsFaces[pointNum])
//         m1, m2, m3 := (n-3)/n, 1.0/n, 2.0/n
//         oldCoords := inputPoints[pointNum]
//         p1 := mulPoint(oldCoords, m1)
//         afp := avgFacePoints[pointNum]
//         p2 := mulPoint(afp, m2)
//         ame := avgMidEdges[pointNum]
//         p3 := mulPoint(ame, m3)
//         p4 := sumPoint(p1, p2)
//         newPoints[pointNum] = sumPoint(p4, p3)
//     }
//     return newPoints
// }

// func switchNums(pointNums [2]int) [2]int {
//     if pointNums[0] < pointNums[1] {
//         return pointNums
//     }
//     return [2]int{pointNums[1], pointNums[0]}
// }

// func cmcSubdiv(inputPoints []Point, inputFaces []Face) ([]Point, []Face) {
//     facePoints := getFacePoints(inputPoints, inputFaces)
//     edgesFaces := getEdgesFaces(inputPoints, inputFaces)
//     edgePoints := getEdgePoints(inputPoints, edgesFaces, facePoints)
//     avgFacePoints := getAvgFacePoints(inputPoints, inputFaces, facePoints)
//     avgMidEdges := getAvgMidEdges(inputPoints, edgesFaces)
//     pointsFaces := getPointsFaces(inputPoints, inputFaces)
//     newPoints := getNewPoints(inputPoints, pointsFaces, avgFacePoints, avgMidEdges)
//     var facePointNums []int
//     nextPointNum := len(newPoints)
//     for _, facePoint := range facePoints {
//         newPoints = append(newPoints, facePoint)
//         facePointNums = append(facePointNums, nextPointNum)
//         nextPointNum++
//     }
//     edgePointNums := make(map[[2]int]int)
//     for edgeNum := range edgesFaces {
//         point1 := edgesFaces[edgeNum].point1
//         point2 := edgesFaces[edgeNum].point2
//         edgePoint := edgePoints[edgeNum]
//         newPoints = append(newPoints, edgePoint)
//         edgePointNums[[2]int{point1, point2}] = nextPointNum
//         nextPointNum++
//     }
//     var newFaces []Face
//     for oldFaceNum, oldFace := range inputFaces {
//         if len(oldFace) == 4 {
//             a, b, c, d := oldFace[0], oldFace[1], oldFace[2], oldFace[3]
//             facePointAbcd := facePointNums[oldFaceNum]
//             edgePointAb := edgePointNums[switchNums([2]int{a, b})]
//             edgePointDa := edgePointNums[switchNums([2]int{d, a})]
//             edgePointBc := edgePointNums[switchNums([2]int{b, c})]
//             edgePointCd := edgePointNums[switchNums([2]int{c, d})]
//             newFaces = append(newFaces, Face{a, edgePointAb, facePointAbcd, edgePointDa})
//             newFaces = append(newFaces, Face{b, edgePointBc, facePointAbcd, edgePointAb})
//             newFaces = append(newFaces, Face{c, edgePointCd, facePointAbcd, edgePointBc})
//             newFaces = append(newFaces, Face{d, edgePointDa, facePointAbcd, edgePointCd})
//         }
//     }
//     return newPoints, newFaces
// }

// func main() {
//     inputPoints := []Point{
//         {-1.0, 1.0, 1.0},
//         {-1.0, -1.0, 1.0},
//         {1.0, -1.0, 1.0},
//         {1.0, 1.0, 1.0},
//         {1.0, -1.0, -1.0},
//         {1.0, 1.0, -1.0},
//         {-1.0, -1.0, -1.0},
//         {-1.0, 1.0, -1.0},
//     }

//     inputFaces := []Face{
//         {0, 1, 2, 3},
//         {3, 2, 4, 5},
//         {5, 4, 6, 7},
//         {7, 0, 3, 5},
//         {7, 6, 1, 0},
//         {6, 1, 2, 4},
//     }

//     outputPoints := make([]Point, len(inputPoints))
//     outputFaces := make([]Face, len(inputFaces))
//     copy(outputPoints, inputPoints)
//     copy(outputFaces, inputFaces)
//     iterations := 1
//     for i := 0; i < iterations; i++ {
//         outputPoints, outputFaces = cmcSubdiv(outputPoints, outputFaces)
//     }
//     for _, p := range outputPoints {
//         fmt.Printf("% .4f\n", p)
//     }
//     fmt.Println()
//     for _, f := range outputFaces {
//         fmt.Printf("%2d\n", f)
//     }
// }
// ^-- This is the go version

// v-- Here is the zig version

const std = @import("std");
const math = std.math;
const ArrayList = std.ArrayList;

const Point = @Vector(3, f32);
const Face = []const u32;
const EdgesFace = struct {
    point1: u32,
    point2: u32,
    face1: u32,
    face2: u32,
    centerPoint: Point,
};
const PointEx = struct {
    p: Point,
    n: u32,
};

fn getFacePoints(allocator: std.mem.Allocator, inputPoints: []Point, inputFaces: []Face) ![]Point {
    var facePoints = try ArrayList(Point).initCapacity(allocator, inputFaces.len);
    for (inputFaces) |face| {
        var facePoint = Point{ 0, 0, 0 };
        for (face) |pointNum| {
            facePoint += inputPoints[pointNum];
        }
        facePoint /= @splat(3, @intToFloat(f32, face.len));
        try facePoints.append(facePoint);
    }
    return facePoints.toOwnedSlice();
}

fn centerPoint(p1: Point, p2: Point) Point {
    return (p1 + p2) / @splat(3, @floatCast(f32, 2));
}

fn getEdgesFaces(allocator: std.mem.Allocator, inputPoints: []Point, inputFaces: []Face) ![]EdgesFace {
    var edges = try ArrayList([3]u32).initCapacity(allocator, inputFaces.len * 4);
    for (inputFaces, 0..) |face, faceNum| {
        const numPoints = face.len;
        for (face, 0..) |pointNum, pointIndex| {
            var point1 = pointNum;
            var point2: u32 = 0;
            if (pointIndex < numPoints - 1) {
                point2 = face[pointIndex + 1];
            } else {
                point2 = face[0];
            }
            if (point1 > point2) {
                var swap = point1;
                point1 = point2;
                point2 = swap;
            }
            try edges.append([3]u32{ point1, point2, @intCast(u32, faceNum) });
        }
    }
    std.sort.sort([3]u32, edges.items, {}, struct {
        fn sort(context: void, a: [3]u32, b: [3]u32) bool {
            _ = context;
            if (a[0] == b[0]) {
                if (a[1] == b[1]) {
                    return a[2] < b[2];
                }
                return a[1] < b[1];
            }
            return a[0] < b[0];
        }
    }.sort);
    var numEdges = edges.items.len;
    var eIndex: usize = 0;
    var mergedEdges = try ArrayList([4]u32).initCapacity(allocator, numEdges);
    while (eIndex < numEdges) : (eIndex += 1) {
        var e1 = edges.items[eIndex];
        if (eIndex < numEdges - 1) {
            var e2 = edges.items[eIndex + 1];
            if (e1[0] == e2[0] and e1[1] == e2[1]) {
                try mergedEdges.append([4]u32{ e1[0], e1[1], e1[2], e2[2] });
                eIndex += 1;
            } else {
                try mergedEdges.append([4]u32{ e1[0], e1[1], e1[2], std.math.maxInt(u32) });
            }
        } else {
            try mergedEdges.append([4]u32{ e1[0], e1[1], e1[2], std.math.maxInt(u32) });
        }
    }
    var edgesCenters = try ArrayList(EdgesFace).initCapacity(allocator, mergedEdges.items.len);
    for (mergedEdges.items) |me| {
        var p1 = inputPoints[me[0]];
        var p2 = inputPoints[me[1]];
        try edgesCenters.append(EdgesFace{
            .point1 = me[0],
            .point2 = me[1],
            .face1 = me[2],
            .face2 = me[3],
            .centerPoint = centerPoint(p1, p2),
        });
    }
    return edgesCenters.items;
}

fn getEdgePoints(allocator: std.mem.Allocator, edgesFaces: []EdgesFace, facePoints: []Point) ![]Point {
    var edgePoints = try ArrayList(Point).initCapacity(allocator, edgesFaces.len);
    for (edgesFaces) |edge| {
        var cp = edge.centerPoint;
        var fp1 = facePoints[edge.face1];
        var fp2: Point = undefined;
        if (edge.face2 == std.math.maxInt(u32)) {
            fp2 = fp1;
        } else {
            fp2 = facePoints[edge.face2];
        }
        var cfp = centerPoint(fp1, fp2);
        try edgePoints.append(centerPoint(cp, cfp));
    }
    return edgePoints.toOwnedSlice();
}

// func getAvgFacePoints(inputPoints []Point, inputFaces []Face, facePoints []Point) []Point {
//     numPoints := len(inputPoints)
//     tempPoints := make([]PointEx, numPoints)
//     for faceNum := range inputFaces {
//         fp := facePoints[faceNum]
//         for _, pointNum := range inputFaces[faceNum] {
//             tp := tempPoints[pointNum].p
//             tempPoints[pointNum].p = sumPoint(tp, fp)
//             tempPoints[pointNum].n++
//         }
//     }
//     avgFacePoints := make([]Point, numPoints)
//     for i, tp := range tempPoints {
//         avgFacePoints[i] = divPoint(tp.p, float64(tp.n))
//     }
//     return avgFacePoints
// }

fn getAvgFacePoints(allocator: std.mem.Allocator, inputPoints: []Point, inputFaces: []Face, facePoints: []Point) ![]Point {
    var tempPoints = try ArrayList(PointEx).initCapacity(allocator, inputPoints.len);
    for (inputPoints) |_| {
        try tempPoints.append(PointEx{ .p = Point{ 0, 0, 0 }, .n = 0 });
    }
    for (inputFaces, 0..) |face, faceNum| {
        var fp = facePoints[faceNum];
        for (face) |pointNum| {
            var tp = tempPoints.items[pointNum].p;
            tempPoints.items[pointNum].p = tp + fp;
            tempPoints.items[pointNum].n += 1;
        }
    }
    var avgFacePoints = try ArrayList(Point).initCapacity(allocator, tempPoints.items.len);
    for (tempPoints.items) |tp| {
        try avgFacePoints.append(tp.p / @splat(3, @intToFloat(f32, tp.n)));
    }
    return avgFacePoints.toOwnedSlice();
}

// func getAvgMidEdges(inputPoints []Point, edgesFaces []Edge) []Point {
//     numPoints := len(inputPoints)
//     tempPoints := make([]PointEx, numPoints)
//     for _, edge := range edgesFaces {
//         centerPoint := edge.centerPoint
//         for _, pointNum := range []int{edge.point1, edge.point2} {
//             tp := tempPoints[pointNum].p
//             tempPoints[pointNum].p = sumPoint(tp, centerPoint)
//             tempPoints[pointNum].n++
//         }
//     }
//     avgMidEdges := make([]Point, len(tempPoints))
//     for i, tp := range tempPoints {
//         avgMidEdges[i] = divPoint(tp.p, float64(tp.n))
//     }
//     return avgMidEdges
// }

fn getAvgMidEdges(allocator: std.mem.Allocator, inputPoints: []Point, edgesFaces: []EdgesFace) ![]Point {
    var tempPoints = try ArrayList(PointEx).initCapacity(allocator, inputPoints.len);
    for (inputPoints) |_| {
        try tempPoints.append(PointEx{ .p = Point{ 0, 0, 0 }, .n = 0 });
    }
    for (edgesFaces) |edge| {
        var edges: []const u32 = &.{ edge.point1, edge.point2 };
        for (edges) |pointNum| {
            var tp = tempPoints.items[pointNum].p;
            tempPoints.items[pointNum].p = tp + edge.centerPoint;
            tempPoints.items[pointNum].n += 1;
        }
    }
    var avgMidEdges = try ArrayList(Point).initCapacity(allocator, tempPoints.items.len);
    for (tempPoints.items) |tp| {
        try avgMidEdges.append(tp.p / @splat(3, @intToFloat(f32, tp.n)));
    }
    return avgMidEdges.toOwnedSlice();
}

// func getPointsFaces(inputPoints []Point, inputFaces []Face) []int {
//     numPoints := len(inputPoints)
//     pointsFaces := make([]int, numPoints)
//     for faceNum := range inputFaces {
//         for _, pointNum := range inputFaces[faceNum] {
//             pointsFaces[pointNum]++
//         }
//     }
//     return pointsFaces
// }

fn getPointsFaces(allocator: std.mem.Allocator, inputPoints: []Point, inputFaces: []Face) ![]u32 {
    var pointsFaces = try ArrayList(u32).initCapacity(allocator, inputPoints.len);
    for (inputPoints) |_| {
        try pointsFaces.append(0);
    }
    for (inputFaces) |face| {
        for (face) |pointNum| {
            pointsFaces.items[pointNum] += 1;
        }
    }
    return pointsFaces.toOwnedSlice();
}

// func getNewPoints(inputPoints []Point, pointsFaces []int, avgFacePoints, avgMidEdges []Point) []Point {
//     newPoints := make([]Point, len(inputPoints))
//     for pointNum := range inputPoints {
//         n := float64(pointsFaces[pointNum])
//         m1, m2, m3 := (n-3)/n, 1.0/n, 2.0/n
//         oldCoords := inputPoints[pointNum]
//         p1 := mulPoint(oldCoords, m1)
//         afp := avgFacePoints[pointNum]
//         p2 := mulPoint(afp, m2)
//         ame := avgMidEdges[pointNum]
//         p3 := mulPoint(ame, m3)
//         p4 := sumPoint(p1, p2)
//         newPoints[pointNum] = sumPoint(p4, p3)
//     }
//     return newPoints
// }

fn getNewPoints(allocator: std.mem.Allocator, inputPoints: []Point, pointsFaces: []u32, avgFacePoints: []Point, avgMidEdges: []Point) ![]Point {
    var newPoints = try ArrayList(Point).initCapacity(allocator, inputPoints.len);
    for (inputPoints, 0..) |point, pointNum| {
        var n = @intToFloat(f32, pointsFaces[pointNum]);
        var m1 = (n - 3) / n;
        var m2 = 1.0 / n;
        var m3 = 2.0 / n;
        var p1 = point * @splat(3, m1);
        var afp = avgFacePoints[pointNum];
        var p2 = afp * @splat(3, m2);
        var ame = avgMidEdges[pointNum];
        var p3 = ame * @splat(3, m3);
        var p4 = p1 + p2;
        try newPoints.append(p4 + p3);
    }
    return newPoints.toOwnedSlice();
}

// func switchNums(pointNums [2]int) [2]int {
//     if pointNums[0] < pointNums[1] {
//         return pointNums
//     }
//     return [2]int{pointNums[1], pointNums[0]}
// }

fn switchNums(pointNums: [2]u32) [2]u32 {
    if (pointNums[0] < pointNums[1]) {
        return pointNums;
    }
    return [2]u32{ pointNums[1], pointNums[0] };
}

// func cmcSubdiv(inputPoints []Point, inputFaces []Face) ([]Point, []Face) {
//     facePoints := getFacePoints(inputPoints, inputFaces)
//     edgesFaces := getEdgesFaces(inputPoints, inputFaces)
//     edgePoints := getEdgePoints(inputPoints, edgesFaces, facePoints)
//     avgFacePoints := getAvgFacePoints(inputPoints, inputFaces, facePoints)
//     avgMidEdges := getAvgMidEdges(inputPoints, edgesFaces)
//     pointsFaces := getPointsFaces(inputPoints, inputFaces)
//     newPoints := getNewPoints(inputPoints, pointsFaces, avgFacePoints, avgMidEdges)
//     var facePointNums []int
//     nextPointNum := len(newPoints)
//     for _, facePoint := range facePoints {
//         newPoints = append(newPoints, facePoint)
//         facePointNums = append(facePointNums, nextPointNum)
//         nextPointNum++
//     }
//     edgePointNums := make(map[[2]int]int)
//     for edgeNum := range edgesFaces {
//         point1 := edgesFaces[edgeNum].point1
//         point2 := edgesFaces[edgeNum].point2
//         edgePoint := edgePoints[edgeNum]
//         newPoints = append(newPoints, edgePoint)
//         edgePointNums[[2]int{point1, point2}] = nextPointNum
//         nextPointNum++
//     }
//     var newFaces []Face
//     for oldFaceNum, oldFace := range inputFaces {
//         if len(oldFace) == 4 {
//             a, b, c, d := oldFace[0], oldFace[1], oldFace[2], oldFace[3]
//             facePointAbcd := facePointNums[oldFaceNum]
//             edgePointAb := edgePointNums[switchNums([2]int{a, b})]
//             edgePointDa := edgePointNums[switchNums([2]int{d, a})]
//             edgePointBc := edgePointNums[switchNums([2]int{b, c})]
//             edgePointCd := edgePointNums[switchNums([2]int{c, d})]
//             newFaces = append(newFaces, Face{a, edgePointAb, facePointAbcd, edgePointDa})
//             newFaces = append(newFaces, Face{b, edgePointBc, facePointAbcd, edgePointAb})
//             newFaces = append(newFaces, Face{c, edgePointCd, facePointAbcd, edgePointBc})
//             newFaces = append(newFaces, Face{d, edgePointDa, facePointAbcd, edgePointCd})
//         }
//     }
//     return newPoints, newFaces
// }

fn cmcSubdiv(allocator: std.mem.Allocator, inputPoints: []Point, inputFaces: []Face) !struct { points: []Point, faces: []Face } {
    var facePoints = try getFacePoints(allocator, inputPoints, inputFaces);
    var edgesFaces = try getEdgesFaces(allocator, inputPoints, inputFaces);
    var edgePoints = try getEdgePoints(allocator, edgesFaces, facePoints);
    var avgFacePoints = try getAvgFacePoints(allocator, inputPoints, inputFaces, facePoints);
    var avgMidEdges = try getAvgMidEdges(allocator, inputPoints, edgesFaces);
    var pointsFaces = try getPointsFaces(allocator, inputPoints, inputFaces);
    var initialNewPoints = try getNewPoints(allocator, inputPoints, pointsFaces, avgFacePoints, avgMidEdges);
    var facePointNums = try ArrayList(u32).initCapacity(allocator, facePoints.len);
    var newPoints = try ArrayList(Point).initCapacity(allocator, initialNewPoints.len);
    try newPoints.appendSlice(initialNewPoints);
    var nextPointNum = newPoints.items.len;
    for (facePoints) |facePoint| {
        try newPoints.append(facePoint);
        try facePointNums.append(@intCast(u32, nextPointNum));
        nextPointNum += 1;
    }
    var edgePointNums = std.AutoHashMap([2]u32, u32).init(allocator);
    for (edgesFaces, 0..) |edgeFace, edgeNum| {
        var point1 = edgeFace.point1;
        var point2 = edgeFace.point2;
        var edgePoint = edgePoints[edgeNum];
        try newPoints.append(edgePoint);
        try edgePointNums.put(switchNums([2]u32{ point1, point2 }), @intCast(u32, nextPointNum));
        nextPointNum += 1;
    }
    var newFaces = try ArrayList(Face).initCapacity(allocator, inputFaces.len);
    for (inputFaces, 0..) |oldFace, oldFaceNum| {
        if (oldFace.len == 4) {
            var a = oldFace[0];
            var b = oldFace[1];
            var c = oldFace[2];
            var d = oldFace[3];
            var facePointAbcd = facePointNums.items[oldFaceNum];
            var edgePointAb = edgePointNums.get(switchNums([2]u32{ a, b })).?;
            var edgePointDa = edgePointNums.get(switchNums([2]u32{ d, a })).?;
            var edgePointBc = edgePointNums.get(switchNums([2]u32{ b, c })).?;
            var edgePointCd = edgePointNums.get(switchNums([2]u32{ c, d })).?;
            try newFaces.append(&.{ a, edgePointAb, facePointAbcd, edgePointDa });
            try newFaces.append(&.{ b, edgePointBc, facePointAbcd, edgePointAb });
            try newFaces.append(&.{ c, edgePointCd, facePointAbcd, edgePointBc });
            try newFaces.append(&.{ d, edgePointDa, facePointAbcd, edgePointCd });
        }
    }
    return .{ .points = try newPoints.toOwnedSlice(), .faces = try newFaces.toOwnedSlice() };
}

test "getFacePoints" {
    var allocator = std.heap.page_allocator;
    var points = [_]Point{
        Point{ -1.0, 1.0, 1.0 },
        Point{ -1.0, -1.0, 1.0 },
        Point{ 1.0, -1.0, 1.0 },
        Point{ 1.0, 1.0, 1.0 },
        Point{ -1.0, 1.0, -1.0 },
        Point{ -1.0, -1.0, -1.0 },
    };
    var constFaces = [_][]const u32{
        &.{ 0, 1, 2, 3 },
        &.{ 0, 1, 5, 4 },
    };
    var faces = try ArrayList(Face).initCapacity(allocator, constFaces.len);
    for (constFaces) |constFace| {
        var face = try ArrayList(u32).initCapacity(allocator, constFace.len);
        for (constFace) |pointNum| {
            try face.append(pointNum);
        }
        try faces.append(face.items);
    }
    var result = try getFacePoints(
        allocator,
        &points,
        faces.items,
    );

    var expected = [_]Point{
        Point{ 0.0, 0.0, 1.0 },
        Point{ -1.0, 0.0, 0.0 },
    };

    try std.testing.expectEqual(expected.len, result.len);
    for (expected, 0..) |expectedPoint, i| {
        try std.testing.expectEqual(expectedPoint, result[i]);
    }
}

test "getEdgesFaces" {
    var allocator = std.heap.page_allocator;
    var points = [_]Point{
        Point{ -1.0, 1.0, 1.0 },
        Point{ -1.0, -1.0, 1.0 },
        Point{ 1.0, -1.0, 1.0 },
        Point{ 1.0, 1.0, 1.0 },
        Point{ -1.0, 1.0, -1.0 },
        Point{ -1.0, -1.0, -1.0 },
    };
    var constFaces = [_][]const u32{
        &.{ 0, 1, 2, 3 },
        &.{ 0, 1, 5, 4 },
    };
    var faces = try ArrayList(Face).initCapacity(allocator, constFaces.len);
    for (constFaces) |constFace| {
        var face = try ArrayList(u32).initCapacity(allocator, constFace.len);
        for (constFace) |pointNum| {
            try face.append(pointNum);
        }
        try faces.append(face.items);
    }
    var result = try getEdgesFaces(
        allocator,
        &points,
        faces.items,
    );

    try std.testing.expectEqual(EdgesFace{ .point1 = 0, .point2 = 1, .face1 = 0, .face2 = 1, .centerPoint = .{ -1.0, 0.0, 1.0 } }, result[0]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 0, .point2 = 3, .face1 = 0, .face2 = std.math.maxInt(u32), .centerPoint = .{ 0.0, 1.0, 1.0 } }, result[1]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 0, .point2 = 4, .face1 = 1, .face2 = std.math.maxInt(u32), .centerPoint = .{ -1.0, 1.0, 0.0 } }, result[2]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 1, .point2 = 2, .face1 = 0, .face2 = std.math.maxInt(u32), .centerPoint = .{ 0.0, -1.0, 1.0 } }, result[3]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 1, .point2 = 5, .face1 = 1, .face2 = std.math.maxInt(u32), .centerPoint = .{ -1.0, -1.0, 0.0 } }, result[4]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 2, .point2 = 3, .face1 = 0, .face2 = std.math.maxInt(u32), .centerPoint = .{ 1.0, 0.0, 1.0 } }, result[5]);
    try std.testing.expectEqual(EdgesFace{ .point1 = 4, .point2 = 5, .face1 = 1, .face2 = std.math.maxInt(u32), .centerPoint = .{ -1.0, 0.0, -1.0 } }, result[6]);
}

test "getPointsFaces" {
    var allocator = std.heap.page_allocator;
    var points = [_]Point{
        Point{ -1.0, 1.0, 1.0 },
        Point{ -1.0, -1.0, 1.0 },
        Point{ 1.0, -1.0, 1.0 },
        Point{ 1.0, 1.0, 1.0 },
        Point{ -1.0, 1.0, -1.0 },
        Point{ -1.0, -1.0, -1.0 },
    };
    var constFaces = [_][]const u32{
        &.{ 0, 1, 2, 3 },
        &.{ 0, 1, 5, 4 },
    };
    var faces = try ArrayList(Face).initCapacity(allocator, constFaces.len);
    for (constFaces) |constFace| {
        var face = try ArrayList(u32).initCapacity(allocator, constFace.len);
        for (constFace) |pointNum| {
            try face.append(pointNum);
        }
        try faces.append(face.items);
    }
    var result = try getPointsFaces(
        allocator,
        &points,
        faces.items,
    );

    _ = result;
}

test "getNewPoints" {
    var allocator = std.heap.page_allocator;
    var points = [_]Point{
        Point{ -1.0, 1.0, 1.0 },
        Point{ -1.0, -1.0, 1.0 },
        Point{ 1.0, -1.0, 1.0 },
        Point{ 1.0, 1.0, 1.0 },
        Point{ -1.0, 1.0, -1.0 },
        Point{ -1.0, -1.0, -1.0 },
    };
    var constFaces = [_][]const u32{
        &.{ 0, 1, 2, 3 },
        &.{ 0, 1, 5, 4 },
    };
    var faces = try ArrayList(Face).initCapacity(allocator, constFaces.len);
    for (constFaces) |constFace| {
        var face = try ArrayList(u32).initCapacity(allocator, constFace.len);
        for (constFace) |pointNum| {
            try face.append(pointNum);
        }
        try faces.append(face.items);
    }
    var pointFaces = try getPointsFaces(
        allocator,
        &points,
        faces.items,
    );
    var facePoints = try getFacePoints(
        allocator,
        &points,
        faces.items,
    );
    var edgesFaces = try getEdgesFaces(
        allocator,
        &points,
        faces.items,
    );
    var avgFacePoints = try getAvgFacePoints(
        allocator,
        &points,
        faces.items,
        facePoints,
    );
    var avgMidEdges = try getAvgMidEdges(
        allocator,
        &points,
        edgesFaces,
    );
    var result = try getNewPoints(
        allocator,
        &points,
        pointFaces,
        avgFacePoints,
        avgMidEdges,
    );

    _ = result;
}

test "cmcSubdiv" {
    var allocator = std.heap.page_allocator;
    var points = [_]Point{
        Point{ -1.0, 1.0, 1.0 },
        Point{ -1.0, -1.0, 1.0 },
        Point{ 1.0, -1.0, 1.0 },
        Point{ 1.0, 1.0, 1.0 },
        Point{ -1.0, 1.0, -1.0 },
        Point{ -1.0, -1.0, -1.0 },
    };
    var constFaces = [_][]const u32{
        &.{ 0, 1, 2, 3 },
        &.{ 0, 1, 5, 4 },
    };
    var faces = try ArrayList(Face).initCapacity(allocator, constFaces.len);
    for (constFaces) |constFace| {
        var face = try ArrayList(u32).initCapacity(allocator, constFace.len);
        for (constFace) |pointNum| {
            try face.append(pointNum);
        }
        try faces.append(face.items);
    }
    var result = try cmcSubdiv(
        allocator,
        &points,
        faces.items,
    );

    _ = result;
}
