const std = @import("std");

const pi: f64 = 3.141592653589793;
const solar_mass: f64 = 4.0 * pi * pi;
const year: f64 = 365.24;
const n_bodies: usize = 5;

const Body = struct {
    x: [3]f64,
    v: [3]f64,
    mass: f64,
};

const bodies = [_]Body{
    // Sun
    Body{
        .x = .{ 0.0, 0.0, 0.0 },
        .v = .{ 0.0, 0.0, 0.0 },
        .mass = solar_mass,
    },
    // Jupiter
    Body{
        .x = .{
            4.84143144246472090e+00,
            -1.16032004402742839e+00,
            -1.03622044471123109e-01,
        },
        .v = .{
            1.66007664274403694e-03 * year,
            7.69901118419740425e-03 * year,
            -6.90460016972063023e-05 * year,
        },
        .mass = 9.54791938424326609e-04 * solar_mass,
    },
    // Saturn
    Body{
        .x = .{
            8.34336671824457987e+00,
            4.12479856412430479e+00,
            -4.03523417114321381e-01,
        },
        .v = .{
            -2.76742510726862411e-03 * year,
            4.99852801234917238e-03 * year,
            2.30417297573763929e-05 * year,
        },
        .mass = 2.85885980666130812e-04 * solar_mass,
    },
    // Uranus
    Body{
        .x = .{
            1.28943695621391310e+01,
            -1.51111514016986312e+01,
            -2.23307578892655734e-01,
        },
        .v = .{
            2.96460137564761618e-03 * year,
            2.37847173959480950e-03 * year,
            -2.96589568540237556e-05 * year,
        },
        .mass = 4.36624404335156298e-05 * solar_mass,
    },
    // Neptune
    Body{
        .x = .{
            1.53796971148509165e+01,
            -2.59193146099879641e+01,
            1.79258772950371181e-01,
        },
        .v = .{
            2.68067772490389322e-03 * year,
            1.62824170038242295e-03 * year,
            -9.51592254519715870e-05 * year,
        },
        .mass = 5.15138902046611451e-05 * solar_mass,
    },
};

pub fn energy(bs: []Body) f64 {
    var e: f64 = 0;
    var delta: [3]f64 = undefined;

    var i: u8 = 0;
    while (i < bs.len) : (i += 1) {
        e += 0.5 * bs[i].mass * (bs[i].v[0] * bs[i].v[0] +
            bs[i].v[1] * bs[i].v[1] +
            bs[i].v[2] * bs[i].v[2]);

        var j: u8 = i + 1;
        while (j < bs.len) : (j += 1) {
            delta[0] = bs[i].x[0] - bs[j].x[0];
            delta[1] = bs[i].x[1] - bs[j].x[1];
            delta[2] = bs[i].x[2] - bs[j].x[2];

            const d_squared = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
            const distance = std.math.sqrt(d_squared);

            e -= (bs[i].mass * bs[j].mass) / distance;
        }
    }
    return e;
}

pub fn offsetMomentum(bs: []Body) void {
    for (bs) |body| {
        for (body.x) |_, k| {
            bs[0].v[k] -= body.v[k] * body.mass / solar_mass;
        }
    }
}

pub fn advance(bs: []Body, dt: f64) void {
    var delta: [3]f64 = undefined;
    var i: u8 = 0;
    while (i < bs.len) : (i += 1) {
        var j: u8 = i + 1;
        while (j < bs.len) : (j += 1) {
            delta[0] = bs[i].x[0] - bs[j].x[0];
            delta[1] = bs[i].x[1] - bs[j].x[1];
            delta[2] = bs[i].x[2] - bs[j].x[2];

            const d_squared = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
            const distance = std.math.sqrt(d_squared);
            const mag = dt / (d_squared * distance);

            bs[i].v[0] -= delta[0] * bs[j].mass * mag;
            bs[i].v[1] -= delta[1] * bs[j].mass * mag;
            bs[i].v[2] -= delta[2] * bs[j].mass * mag;

            bs[j].v[0] += delta[0] * bs[i].mass * mag;
            bs[j].v[1] += delta[1] * bs[i].mass * mag;
            bs[j].v[2] += delta[2] * bs[i].mass * mag;
        }
    }

    for (bs) |*body| {
        body.x[0] += dt * body.v[0];
        body.x[1] += dt * body.v[1];
        body.x[2] += dt * body.v[2];
    }
}

pub fn main() void {
    var bs = bodies;

    offsetMomentum(bs[0..]);
    std.debug.warn("{d:.9}\n", .{energy(bs[0..])});

    var n: usize = 0;
    while (n < 50000000) : (n += 1) {
        advance(bs[0..], 0.01);
    }

    std.debug.warn("{d:.9}\n", .{energy(bs[0..])});
}
