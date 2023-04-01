const std = @import("std");

fn ReactiveModel(comptime T: type) type {
    return struct {
        const Self = @This();
        value: T,
        fn set(self: *Self, value: T) void {
            self.value = value;
        }
        fn init(value: T) Self {
            return Self{ .value = value };
        }
        fn register(self: Self, comptime system: type) void {
            std.debug.print("registering {*}\n", .{ @typeName(system) });
            const Args = @typeInfo(@TypeOf(system.listen)).Fn.params[0].type.?;
            const args = Args{ .foo = self.value.foo };
            system.listen(args);
        }        
    };
}

const test_model = struct { foo: u32 };

const test_system = struct {
    foo: u32,
    fn listen(thinger: @This()) void {
        std.debug.print("foo = {}\n", .{thinger.foo});
    }
};

test "reactive model" {
    var model = ReactiveModel(test_model).init(.{ .foo = 0 });
    model.set(.{ .foo = 1 });
    try std.testing.expect(model.value.foo == 1);
}

test "register" {
    var model = ReactiveModel(test_model).init(.{ .foo = 0 });
    model.register(test_system);
}