# Ninja Physics

Ninja Physics is a physics plugin for LesserPanda engine, and is inspired by [articles](http://www.metanetsoftware.com/technique/tutorialA.html) of Metanet Software.

## Features

- Spatial partitioning based broad phase
- Better SAT implementation (Uses modified [SAT.js](https://github.com/jriecken/sat-js))
- Tile collision (coming soon)

## Install

Simply download "ninja-physics.js" and put it inside "game/plugins" folder. Don't forget
to "require" it anywhere you want.

This plugin ONLY changes the implementation and add more features, the API is almost the same of
official physics module.

## Usage

### Basic

You can use `Body` and `World` the same as before.

```javascript
// Create a body like you've always done
var circleBody = new game.Body({
  position: { x: 120, y: 120 },
  shape: new game.Circle(40),
  mass: 1, // will fall
  collisionGroup: 1,
  collideAgainst: [2]
});

var rectBody = new game.Body({
  position: { x: 100, y: 100 },
  shape: new game.Rectangle(40, 40),
  collisionGroup: 2
});

// Add them to world and it just works
circleBody.addTo(game.scene.world);
rectBody.addTo(game.scene.world);
```

### Rotate Body

Body is not limited to use AABB for collision detection, you can
rotate it and still get the right behavior.

```javascript
// Create a body
var body = new game.Body({
  position: { x: 100, y: 100 },
  shape: new game.Rectangle(40, 40)
});

// Rotate it a little bit (45 degrees)
body.shape.setAngle(Math.PI * 0.25);
```

### Manually Check Body Overlapping

There're some utility functions you can use to detect
body overlapping, which will give you a detailed result.

```javascript
// Create 2 circles
var circle1 = new game.Body({
    position: { x: 0, y: 0 },
    shape: new game.Circle(20)
});
var circle2 = new game.Body({
    position: { x: 30, y: 0 },
    shape: new game.Circle(20)
});

// Create a response object to receive collision result
var response = new game.Response();

// Test collision between them
var collided = game.testCircleCircle(circle1, circle2, response);

collided === true; // The two circles are overlapping
response.overlap === 10; // There is 10px overlapping between them
response.overlapV.x === 10; // They are only overlapping in x-axis
response.overlapV.y === 0; // No overlapping in y-axis
```

### Advance

#### Collision group explaination

Collision groups perform a little different from Panda built-in physics
module. Forturnately, the difference is very small.

I'm going to use "Breakout" as the example to explain how you can define
collision groups.

```javascript
// It is always recommended to define groups as constants instead of
// using number directly.
var GROUPS = {
  WALL: 0,   // Wall, solid and no need to customize its `collide` method
  BLOCK: 1,   // Blocks, use `collide` to spawn items and update score
  BALL: 2,    // Ball, will bounce back from WALL and BLOCK
  PADDLE: 3   // Paddle, control how ball is going to move on collision
};

// Wall does not have `collideAgainst`,
// so `collide` method won't even be called once and do not receive
// any responses.
//
// Which means that it is completely SOLID and
// sometimes people called it "FIXED".
var wall = new game.Body({
  collisionGroup: GROUPS.WALL
});

// Blocks will collide with ball.
var block = new game.Body({
  collisionGroup: GROUPS.BLOCK,
  collideAgainst: [GROUPS.BALL],
  collide: function(ball) {
    // Get score, remove self from world
    block.remove();

    // Do not move
    return false;
  }
});

var ball = new game.Body({
  collisionGroup: GROUPS.BALL,
  collideAgainst: [GROUPS.WALL, GROUPS.BLOCK],
  collide: function(other, response) {
    // Bounce back
    var n = response.overlapN.clone();
    this.velocity.subtract(n.multiply(2 * this.velocity.dot(n)));
    // Receive response
    return true;
  }
});

var paddle = new game.Body({
  collisionGroup: GROUPS.PADDLE,
  collideAgainst: [GROUPS.BALL],
  collide: function(ball, response) {
    // Bounce the ball back based on its position
    // ...

    // Do not move
    return false;
  }
});
```

## API Doc

### Polygon

There is a new `Polygon` shape added, and `Rectangle` shape
will be automatically converted to a `Polygon`.

```javascript
/**
 * @constructor
 * points {Array} points Array of vectors representing the original points of the polygon.
 */
var polygon = new game.Polygon(points);

/**
 * Set the current rotation angle (in radians).
 * @method
 */
polygon.setAngle(angle);
```

### Response

You have to create Response objects to receive collision
result from collision testing functions. The `collide` method
of `Body` will also receive a result as second parameter.

```javascript
/**
 * @constructor
 */
var response = new game.Response();

/**
 * The first object in the collision.
 * @property {game.Body}
 */
response.a;

/**
 * The second object in the collison.
 * @property {game.Body}
 */
response.b;

/**
 * Magnitude of the overlap on the shortest colliding axis.
 * @property {Number}
 */
response.overlap;

/**
 * The shortest colliding axis (unit-vector)
 * @property {game.Vector}
 */
response.overlapN;

/**
 * The overlap vector (i.e. overlapN.scale(overlap, overlap)).
 * If this vector is subtracted from the position of a, a and b will no longer be colliding.
 * @property {game.Vector}
 */
response.overlapV;

/**
 * Whether the first object is completely inside the second.
 * @property {Boolean}
 */
response.aInB;

/**
 * Whether the second object is completely inside the first.
 * @property {Boolean}
 */
response.bInA;

/**
 * Clear the response so that it is ready to be reused for another collision test.
 * @method
 */
response.clear();
```

### Collision Testing Methods

```javascript
/**
 * Check if a point is inside a circle.
 * @param {game.Vector} p The point to test
 * @param {game.Circle} c The circle to test
 * @return {Boolean} true if the point is inside the circle, false if it is not
 */
game.pointInCircle(p, c)

/**
 * Check if a point is inside a convex polygon.
 * @param {game.Vector} p The point to test
 * @param {game.Polygon} poly The polygon to test
 * @return {Boolean} true if the point is inside the polygon, false if it is not
 */
game.pointInPolygon(p, poly)

/**
 * Check if two circles collide.
 * @param {game.Body} a The first circle body
 * @param {game.Body} b The second circle body
 * @param {game.Response=} response Response object (optional) that will be populated if
 *   the circles intersect
 * @return {Boolean} true if the circles intersect, false if they don't
 */
game.testCircleCircle(a, b, response)

/**
 * Check if a polygon and a circle collide.
 * @param {game.Polygon} polygon The polygon
 * @param {game.Circle} circle The circle
 * @param {game.Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
game.testPolygonCircle(polygon, circle, response)

/**
 * Check if a circle and a polygon collide.
 *
 * @param {game.Circle} circle The circle
 * @param {game.Polygon} polygon The polygon
 * @param {game.Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
game.testCirclePolygon(circle, polygon, response)

/**
 * Checks whether polygons collide.
 * @param {game.Polygon} a The first polygon
 * @param {game.Polygon} b The second polygon
 * @param {game.Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
game.testPolygonPolygon(a, b, response)
```

## ChangeLog

### v0.1.0

- Polygon shape
- Extend `Vector` class
- Collision solver with more detailed `response` feedback

---

The MIT License (MIT)

Copyright (c) 2015 Sean Bohan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
