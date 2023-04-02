const std = @import("std");

// Given a type T, allow systems to register to changes of the model.
// If a system requests a pointer to a member of the model, we should inform
// all down-stream systems of a possible change.
fn ReactiveModel(comptime T: type, comptime systems: anytype) type {
    return struct {
        const Self = @This();
        model: T,
        fn init(model: T) Self {
            return Self{ .model = model };
        }
        fn set(self: *Self, partial_model: anytype) void {
            inline for (@typeInfo(@TypeOf(partial_model)).Struct.fields) |field| {
                @field(self.model, field.name) = @field(partial_model, field.name);
                inline for (systems) |system| {
                    if (@hasField(system, field.name)) {
                        const SystemValues = @typeInfo(@TypeOf(system.listen)).Fn.params[0].type.?;
                        var system_values: SystemValues = undefined;
                        inline for (@typeInfo(SystemValues).Struct.fields) |system_field| {
                            switch (@typeInfo(system_field.type)) {
                                .Pointer => {
                                    @field(system_values, system_field.name) = &@field(self.model, system_field.name);
                                },
                                else => {
                                    @field(system_values, system_field.name) = @field(self.model, system_field.name);
                                },
                            }
                        }
                        system.listen(system_values);
                    }
                }
            }
        }
    };
}

const test_model = struct {
    foo: u32,
    bar: f32,
};

const test_system = struct {
    foo: u32,
    bar: *f32,
    fn listen(thinger: @This()) void {
        thinger.bar.* = @intToFloat(f32, thinger.foo);
        std.debug.print("bar = {}\n", .{thinger.bar.*});
    }
};

test "reactive model" {
    var model = ReactiveModel(test_model, .{}).init(.{ .foo = 0, .bar = 0.0 });
    model.set(.{ .foo = 1, .bar = 1.0 });
    try std.testing.expect(model.model.foo == 1);
}

test "register" {
    var model = ReactiveModel(test_model, .{ test_system })
        .init(.{ .foo = 123, .bar = 456.0 });
    model.set(.{ .foo = 3, .bar = 1.0 });
    std.debug.print("foo = {}\n", .{model.model.foo});
}
