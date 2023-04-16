const std = @import("std");

// Given a type T, allow systems to register to changes of the model.
// If a system requests a pointer to a member of the model, we should inform
// all down-stream systems of a possible change.
fn ResponsiveModel(comptime systems: anytype) type {
    const T = T: {
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
        break :T @Type(t_type_info);
    };
    return struct {
        const Self = @This();
        state: T,
        fn init(state: T) Self {
            return Self{ .state = state };
        }
        fn retrieveRespondersOnField(comptime args: struct {
            responders: []const type,
            changed_field: std.builtin.Type.StructField,
            ignore: ?type = null,
        }) []const type {
            comptime var next_responders: []const type = args.responders;
            inline for (systems) |system| {
                if (args.ignore) |ignore|
                    if (system == ignore) continue;
                if (has_field: {
                    const SystemValues = @typeInfo(@TypeOf(system.respond)).Fn.params[0].type.?;
                    inline for (@typeInfo(SystemValues).Struct.fields) |next_field| {
                        if (comptime std.mem.eql(u8, next_field.name, args.changed_field.name)) {
                            const next_field_type_info = @typeInfo(next_field.type);
                            if (next_field_type_info != .Pointer or next_field_type_info.Pointer.size != .One)
                                break :has_field true;
                        }
                    } else break :has_field false;
                } and !already_next: {
                    inline for (next_responders) |next_responder| {
                        if (next_responder == system) break :already_next true;
                    }
                    break :already_next false;
                }) {
                    next_responders = next_responders ++ .{system};
                }
            }
            return next_responders;
        }
        fn set(self: *Self, partial_model: anytype) void {
            comptime var responders: []const type = &[_]type{};
            inline for (@typeInfo(@TypeOf(partial_model)).Struct.fields) |field| {
                @field(self.state, field.name) = @field(partial_model, field.name);
                responders = retrieveRespondersOnField(.{
                    .responders = responders,
                    .changed_field = field,
                });
            }
            inline while (responders.len > 0) {
                comptime var next_responders: []const type = &.{};
                inline for (responders) |responder| {
                    const SystemValues = @typeInfo(@TypeOf(responder.respond)).Fn.params[0].type.?;
                    var system_values: SystemValues = undefined;
                    inline for (@typeInfo(SystemValues).Struct.fields) |field| {
                        const field_type_info = @typeInfo(field.type);
                        if (field_type_info != .Pointer or field_type_info.Pointer.size != .One) {
                            @field(system_values, field.name) = @field(self.state, field.name);
                        } else {
                            @field(system_values, field.name) = &@field(self.state, field.name);
                            next_responders = retrieveRespondersOnField(.{
                                .responders = next_responders,
                                .changed_field = field,
                                .ignore = responder,
                            });
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
            source: u21,
            list: *List,
            fn respond(self: @This()) void {
                self.list.append(self.source) catch unreachable;
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
        .source = undefined,
        .list = std.ArrayList(u21).init(std.testing.allocator),
        .count = 0,
    });
    defer model.state.list.deinit();

    model.set(.{ .source = '☔' });
    model.set(.{ .source = '😀' });

    try std.testing.expectEqual(model.state.list.items[0], '☔');
    try std.testing.expectEqual(model.state.list.items[1], '😀');
    try std.testing.expectEqual(model.state.count, 2);
}

test "Detect mutation of slices" {
    const List = std.ArrayList(u21);

    const AddElementToListSystem = struct {
        first_value: ?u21,
        list: *List,
        fn respond(self: @This()) void {
            self.list.items[0] = self.first_value.?;
        }
    };

    const CountListSystem = struct {
        list: List,
        count: *usize,
        fn respond(self: @This()) void {
            self.count.* = self.list.items.len;
        }
    };

    var model = ResponsiveModel(.{
        AddElementToListSystem,
        CountListSystem,
    }).init(undefined);

    var list = List.init(std.testing.allocator);
    defer list.deinit();

    try list.append('☔');
    try list.append('😀');

    model.set(.{ .list = list, .first_value = '🌧' });

    try std.testing.expectEqual(model.state.list.items[0], '🌧');
    try std.testing.expectEqual(model.state.list.items[1], '😀');
    try std.testing.expectEqual(model.state.count, 2);
}
