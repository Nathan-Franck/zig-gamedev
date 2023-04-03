const std = @import("std");

// Given a type T, allow systems to register to changes of the model.
// If a system requests a pointer to a member of the model, we should inform
// all down-stream systems of a possible change.
fn ResponsiveModel(comptime systems: anytype) type {
    const T = blk: {
        comptime var fields: []const std.builtin.Type.StructField = &.{};
        inline for (systems) |system| {
            const SystemValues = @typeInfo(@TypeOf(system.respond)).Fn.params[0].type.?;
            system_fields: inline for (@typeInfo(SystemValues).Struct.fields) |field| {
                const field_type_info = @typeInfo(field.type);
                const field_type = if (field_type_info != .Pointer or field_type_info.Pointer.size != .One)
                    field.type
                else
                    field_type_info.Pointer.child;
                inline for (fields) |existing_field|
                    if (std.mem.eql(u8, existing_field.name, field.name)) {
                        if (existing_field.type != field_type)
                            @compileError("field " ++ field.name ++ " already exists with a different type");
                        continue :system_fields;
                    };
                fields = fields ++ .{.{
                    .default_value = null,
                    .is_comptime = false,
                    .alignment = 0,
                    .name = field.name,
                    .type = field_type,
                }};
            }
        }
        const t_type_info: std.builtin.Type = .{ .Struct = .{
            .fields = fields,
            .decls = &[_]std.builtin.Type.Declaration{},
            .layout = .Auto,
            .is_tuple = false,
        } };
        break :blk @Type(t_type_info);
    };
    return struct {
        const Self = @This();
        state: T,
        // listeners:
        fn init(state: T) Self {
            return Self{ .state = state };
        }
        fn set(self: *Self, partial_model: anytype) void {
            comptime var responders: []const type = &[_]type{};
            inline for (@typeInfo(@TypeOf(partial_model)).Struct.fields) |field| {
                @field(self.state, field.name) = @field(partial_model, field.name);
                inline for (systems) |system| {
                    if (@hasField(system, field.name)) {
                        const already_next = blk: {
                            inline for (responders) |next_responder| {
                                if (next_responder == system) break :blk true;
                            }
                            break :blk false;
                        };
                        if (!already_next) {
                            responders = responders ++ &[_]type{system};
                        }
                    }
                }
            }
            inline while (responders.len > 0) {
                comptime var next_responders: []const type = &.{};
                inline for (responders) |responder| {
                    const SystemValues = @typeInfo(@TypeOf(responder.respond)).Fn.params[0].type.?;
                    var system_values: SystemValues = undefined;
                    inline for (@typeInfo(SystemValues).Struct.fields) |system_field| {
                        const field_type_info = @typeInfo(system_field.type);
                        if (field_type_info != .Pointer or field_type_info.Pointer.size != .One) {
                            @field(system_values, system_field.name) = @field(self.state, system_field.name);
                        } else {
                            @field(system_values, system_field.name) = &@field(self.state, system_field.name);
                            inline for (systems) |system| {
                                if (system == responder) continue;
                                if (@hasField(system, system_field.name)) {
                                    const already_next = blk: {
                                        inline for (next_responders) |next_responder| {
                                            if (next_responder == system) break :blk true;
                                        }
                                        break :blk false;
                                    };
                                    if (!already_next) {
                                        next_responders = next_responders ++ .{system};
                                    }
                                }
                            }
                        }
                    }
                    responder.respond(system_values);
                }
                responders = next_responders;
            }
        }
    };
}

test "State composed from system fields" {
    var model = ResponsiveModel(.{
        struct { // Dummy system declares a field that becomes the model's state.
            initial_value: []const u8,
            fn respond(_: @This()) void {}
        },
    }).init(undefined);

    model.set(.{ .initial_value = "hello" });

    try std.testing.expectEqualStrings(model.state.initial_value, "hello");
}

test "Systems can depend on each other" {
    var model = ResponsiveModel(.{
        struct { // First system takes a value and computes a derived value.
            initial_value: []const u8,
            intermediate_value: *usize,
            fn respond(self: @This()) void {
                self.intermediate_value.* = self.initial_value.len;
            }
        },
        struct { // Second system takes the derived value and computes a final value.
            intermediate_value: usize,
            final_value: *usize,
            fn respond(self: @This()) void {
                self.final_value.* = self.intermediate_value * 2;
            }
        },
    }).init(undefined);

    model.set(.{ .initial_value = "hello" });

    try std.testing.expectEqualDeep(model.state, .{
        .initial_value = "hello",
        .intermediate_value = 5,
        .final_value = 10,
    });
}

test "Systems can compose long pipeline of data" {
    var model = ResponsiveModel(.{
        struct { // First system takes a value and computes a derived value.
            initial_value: []const u8,
            intermediate_value: *usize,
            fn respond(self: @This()) void {
                self.intermediate_value.* = self.initial_value.len;
            }
        },
        struct { // Second system takes the derived value and computes a final value.
            intermediate_value: usize,
            final_value: *usize,
            fn respond(self: @This()) void {
                self.final_value.* = self.intermediate_value * 2;
            }
        },
        struct { // Third system takes the final value and computes a FINAL final value.
            final_value: usize,
            final_final_value: *usize,
            fn respond(self: @This()) void {
                self.final_final_value.* = self.final_value * 2;
            }
        },
    }).init(undefined);

    model.set(.{ .initial_value = "hello" });

    try std.testing.expectEqualDeep(model.state, .{
        .initial_value = "hello",
        .intermediate_value = 5,
        .final_value = 10,
        .final_final_value = 20,
    });
}

test "Array List" {
    var list = std.ArrayList(u21).init(std.testing.allocator);
    defer list.deinit();
    try list.append('☔');

    try std.testing.expect(list.items.len == 1);
}

test "Array list can be used in ResponsiveModel" {
    const List = std.ArrayList(u21);
    var model = ResponsiveModel(.{
        struct { // take value and append to list
            source: *?u21,
            list: *List,
            fn respond(self: @This()) void {
                self.list.append(self.source.*.?) catch unreachable;
                self.source.* = null;
            }
        },
        struct { // count list into new field
            list: List,
            count: *usize,
            fn respond(self: @This()) void {
                self.count.* = self.list.items.len;
            }
        },
    }).init(.{
        .source = null,
        .list = std.ArrayList(u21).init(std.testing.allocator),
        .count = 0,
    });
    defer model.state.list.deinit();

    model.set(.{ .source = '☔' });
    model.set(.{ .source = '😀' });

    try std.testing.expectEqual(model.state.list.items[0], '☔');
    try std.testing.expectEqual(model.state.list.items[1], '😀');
    try std.testing.expectEqual(model.state.source, null);
    try std.testing.expectEqual(model.state.count, 2);
}

test "Detect mutation of slices" {
    const List = std.ArrayList(u21);
    var list = List.init(std.testing.allocator);
    defer list.deinit();

    try list.append('☔');
    try list.append('😀');

    var model = ResponsiveModel(.{
        struct { // change first value of list
            first_value: ?u21,
            list: List,
            fn respond(self: @This()) void {
                self.list.items[0] = self.first_value.?;
            }
        },
        struct { // count list into new field
            list: List,
            count: *usize,
            fn respond(self: @This()) void {
                self.count.* = self.list.items.len;
            }
        },
    }).init(.{
        .first_value = null,
        .list = list,
        .count = 0,
    });

    model.set(.{ .first_value = '🌧' });

    try std.testing.expectEqual(model.state.list.items[0], '🌧');
    try std.testing.expectEqual(model.state.list.items[1], '😀');
    // try std.testing.expectEqual(model.state.count, 2); CURRENTLY THIS DOESN'T WORK
}

// test "Modify the value referenced in a const pointer" {
//     var string: []const u8 = "hello";
//     var ptr: *const u8 = &string[0];
//     ptr.* = 'b';
//     try std.testing.expectEqualStrings(string, "bello");
// }
