const std = @import("std");
const math = std.math;
const ArrayList = std.ArrayList;

const vec3 = @Vector(3, f32);

const meshBlueprint = .{
    .vertices = [_]vec3 {
        .{ -1.0, -1.0, -1.0 },
        .{ 1.0, -1.0, -1.0 },
        .{ 1.0, 1.0, -1.0 },
        .{ -1.0, 1.0, -1.0 },
        .{ -1.0, -1.0, 1.0 },
        .{ 1.0, -1.0, 1.0 },
        .{ 1.0, 1.0, 1.0 },
        .{ -1.0, 1.0, 1.0 },
    },
    .faces = [_][4]u32 {
        .{ 0, 1, 2, 3 },
        .{ 4, 5, 6, 7 },
        .{ 7, 6, 2, 3 },
        .{ 4, 7, 3, 0 },
        .{ 5, 4, 0, 1 },
        .{ 6, 5, 1, 2 },
    },
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var allocator = gpa.allocator();

    var vertices = try ArrayList(vec3).initCapacity(allocator, 8);
    for (meshBlueprint.vertices) |v| {
        try vertices.append(v);
    }
    defer vertices.deinit();

    var faces = try ArrayList(ArrayList(u32)).initCapacity(allocator, 6);

    for (meshBlueprint.faces) |f| {
        var face = try ArrayList(u32).initCapacity(allocator, 4);
        try face.appendSlice(&f);
        try faces.append(face);
    }
    defer faces.deinit();
    defer for (faces.items) |f| {
        f.deinit();
    };

    var subdivided = try subdivideMesh(allocator, vertices.items[0..], faces.items[0..]);
    defer subdivided.deinit();

    std.debug.print("Vertices:\n", .{});
    for (subdivided.vertices.items) |v| {
        std.debug.print("{?}\n", .{v});
    }
    std.debug.print("Faces:\n", .{});
    for (subdivided.faces.items) |f| {
        std.debug.print("{any}\n", .{f.items});
    }
}

const Mesh = struct {
    vertices: ArrayList(vec3),
    faces: ArrayList(ArrayList(usize)),
    fn init(allocator: std.mem.Allocator) @This() {
        return @This(){
            .vertices = ArrayList(vec3).init(allocator),
            .faces = ArrayList(ArrayList(usize)).init(allocator),
        };
    }
    fn deinit(self: *@This()) void {
        self.vertices.deinit();
        self.faces.deinit();
    }
};

fn subdivideMesh(
    allocator: std.mem.Allocator,
    vertices: []const vec3,
    faces: []const ArrayList(u32),
) !Mesh {
    var newMesh = Mesh.init(allocator);

    var edgePoints = std.AutoHashMap(struct { u32, u32 }, usize).init(allocator);
    defer edgePoints.deinit();

    // This is the part that's not working
    // for (faces) |f| {
    //     var facePoint = vec3{ 0.0, 0.0, 0.0 };
    //     for (f.items) |v| {
    //         facePoint += vertices[v];
    //     }
    //     facePoint /= @splat(3, @intToFloat(f32, f.items.len));

    //     try newMesh.vertices.append(facePoint);

    //     var newFace = try ArrayList(usize).initCapacity(allocator, f.items.len);
    //     defer newFace.deinit();

    //     for (f.items, 0..) |v1, i| {
    //         const v2 = f.items[(i + 1) % f.items.len];

    //         const edgeKey = if (v1 < v2) .{ v1, v2 } else .{ v2, v1 };

    //         if (edgePoints.get(edgeKey)) |edge| {
    //             try newFace.append(edge);
    //         } else {
    //             var oppositeFacePoint = oppositeFacePoint: {
    //                 var oppositeFace = faces[getOppositeFaceIndex(f, faces, v1, v2)];
    //                 var x = vec3{ 0.0, 0.0, 0.0 };
    //                 for (oppositeFace.items) |v| {
    //                     x += vertices[v];
    //                 }
    //                 x /= @splat(3, @intToFloat(f32, f.items.len));
    //                 break :oppositeFacePoint x;
    //             };
    //             const me = (vertices[v1] + vertices[v2]) / @splat(3, @as(f32, 2.0));
    //             const ep = (facePoint + oppositeFacePoint + me * @splat(3, @as(f32, 2.0))) / @splat(3, @as(f32, 4.0));

    //             try newMesh.vertices.append(ep);
    //             try edgePoints.put(edgeKey, newMesh.vertices.items.len - 1);

    //             try newFace.append(newMesh.vertices.items.len - 1);
    //         }
    //     }

    // This implementation makes more sense!
    for (faces) |f| {
        var facePoint = vec3{ 0.0, 0.0, 0.0 };
        for (f.items) |v| {
            facePoint += vertices[v];
        }
        facePoint /= @splat(3, @intToFloat(f32, f.items.len));

        try newMesh.vertices.append(facePoint);

        for (f.items, 0..) |v1, i| {
            const v2 = f.items[(i + 1) % f.items.len];

            const edgeKey = if (v1 < v2) .{ v1, v2 } else .{ v2, v1 };

            if (edgePoints.get(edgeKey)) |edge| {
                try newMesh.faces.items[newMesh.faces.items.len - 1].append(edge);
            } else {
                var oppositeFacePoint = oppositeFacePoint: {
                    var oppositeFace = faces[getOppositeFaceIndex(f, faces, v1, v2)];
                    var x = vec3{ 0.0, 0.0, 0.0 };
                    for (oppositeFace.items) |v| {
                        x += vertices[v];
                    }
                    x /= @splat(3, @intToFloat(f32, f.items.len));
                    break :oppositeFacePoint x;
                };
                const me = (vertices[v1] + vertices[v2]) / @splat(3, @as(f32, 2.0));
                const ep = (facePoint + oppositeFacePoint + me * @splat(3, @as(f32, 2.0))) / @splat(3, @as(f32, 4.0));

                try newMesh.vertices.append(ep);
                try edgePoints.put(edgeKey, newMesh.vertices.items.len - 1);

                var newFace = try ArrayList(usize).initCapacity(allocator, 4);

                try newFace.append(v1);
                try newFace.append(newMesh.vertices.items.len - 1);
                try newFace.append(v2);
                try newFace.append(newMesh.vertices.items.len - 1);

                try newMesh.faces.append(newFace);
            }
        }
    }

    for (vertices, 0..) |p, i| {
        var f = vec3{ 0.0, 0.0, 0.0 };
        var r = vec3{ 0.0, 0.0, 0.0 };
        var n: i32 = 0;

        for (faces, 0..) |f1, faceIndex| {
            if (indexOfI: {
                for (f1.items, 0..) |v, f1Index| {
                    if (v == i) break :indexOfI f1Index;
                }
                break :indexOfI null;
            }) |indexOfI| {
                f += newMesh.vertices.items[faceIndex];
                const next = f1.items[(indexOfI + 1) % f1.items.len];
                const prev = f1.items[(indexOfI + f1.items.len - 1) % f1.items.len];
                r += (p + vertices[next]) / @splat(3, @as(f32, 2.0));
                r += (p + vertices[prev]) / @splat(3, @as(f32, 2.0));
                n += 1;
            }
        }

        f /= @splat(3, @intToFloat(f32, n));
        r /= @splat(3, @intToFloat(f32, n * 2));
        const vp = (f + r * @splat(3, @as(f32, 2.0)) + p * (@splat(3, @intToFloat(f32, n) - 3))) / @splat(3, @intToFloat(f32, n));

        try newMesh.vertices.append(vp);
    }

    return newMesh;
}

pub fn getOppositeFaceIndex(
    face: ArrayList(u32),
    faces: []const ArrayList(u32),
    v1: u32,
    v2: u32,
) usize {
    for (faces, 0..) |f, i| {
        if (containsV1: {
            for (f.items) |v| {
                if (v == v1) break :containsV1 true;
            }
            break :containsV1 false;
        } and containsV2: {
            for (f.items) |v| {
                if (v == v2) break :containsV2 true;
            }
            break :containsV2 false;
        } and std.meta.eql(f, face)) {
            return i;
        }
    }
    unreachable;
}

const testing = @import("std").testing;

test "getOppositeFaceIndex" {
}
