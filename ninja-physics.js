/**
 * NinjaPhysics for LesserPanda engine
 * @version 1.0.0
 * @author Sean Bohan
 */
import Vector from 'engine/vector';
import Timer from 'engine/timer';
import Scene from 'engine/scene';
import { Rectangle, Circle } from 'engine/physics';
import { Body as BaseBody } from 'engine/physics';

Object.assign(Vector.prototype, {
  /**
   * Change this vector to be perpendicular to what it was before. (Effectively
   * roatates it 90 degrees in a clockwise direction)
   * @return {Vector} This for chaining
   */
  perp: function perp() {
    var x = this.x;
    this.x = this.y;
    this.y = -x;
    return this;
  },
  /**
   * Reverse this vector
   * @return {Vector} This for chaining
   */
  reverse: function reverse() {
    this.x = -this.x;
    this.y = -this.y;
    return this;
  },
  /**
   * Project this vector on to another vector.
   * @param {Vector} other The vector to project onto
   * @return {Vector} This for chaining
   */
  project: function project(other) {
    var amt = this.dot(other) / other.squaredLength();
    this.x = amt * other.x;
    this.y = amt * other.y;
    return this;
  },
  /**
   * Project this vector onto a vector of unit length. This is slightly more efficient
   * than `project` when dealing with unit vectors.
   * @param {Vector} other The unit vector to project onto
   * @return {Vector} This for chaining
   */
  projectN: function projectN(other) {
    var amt = this.dot(other);
    this.x = amt * other.x;
    this.y = amt * other.y;
    return this;
  },
  /**
   * Reflect this vector on an arbitrary axis.
   * @param {Vector} axis The vector representing the axis
   * @return {Vector} This for chaining
   */
  reflect: function reflect(axis) {
    var x = this.x;
    var y = this.y;
    this.project(axis).multiply(2);
    this.x -= x;
    this.y -= y;
    return this;
  },
  /**
   * Reflect this vector on an arbitrary axis (represented by a unit vector). This is
   * slightly more efficient than `reflect` when dealing with an axis that is a unit vector.
   * @param {Vector} axis The unit vector representing the axis
   * @return {Vector} This for chaining
   */
  reflectN: function reflectN(axis) {
    var x = this.x;
    var y = this.y;
    this.projectN(axis).multiply(2);
    this.x -= x;
    this.y -= y;
    return this;
  },
  /**
   * Returns positive if vector is clockwise of this vector
   * @return {Number} 1: CW, -1: CCW
   */
  sign: function sign(vector) {
    return (this.y * vector.x > this.x * vector.y) ? -1 : 1;
  },
});

Vector.ZERO = new Vector();
Vector.UP = new Vector(0, -1);
Vector.DOWN = new Vector(0, 1);
Vector.LEFT = new Vector(-1, 0);
Vector.RIGHT = new Vector(1, 0);

// TODO: finish polygon related methods
class Polygon {
  constructor(points = []) {
    this.points = [];
    this.calcPoints = [];
    this.edges = [];
    this.normals = [];
    this._rotation = 0;
    this.offset = new Vector();
    this.setPoints(points);
  }

  /**
   * Set the points of the polygon.
   * @param {Array<Vector>=} points An array of vectors representing the points in the polygon,
   *   in clockwise order
   * @return {Polygon} This for chaining
   */
  setPoints(points) {
    // Only re-allocate if this is a new polygon or the number of points has changed.
    var lengthChanged = !this.points || this.points.length !== points.length;
    if (lengthChanged) {
      var calcPoints = this.calcPoints = [];
      var edges = this.edges = [];
      var normals = this.normals = [];
      // Allocate the vector arrays for the calculated properties
      for (var i = 0, len = points.length; i < len; i++) {
        calcPoints.push(new Vector());
        edges.push(new Vector());
        normals.push(new Vector());
      }
    }
    this.points = points;
    this._recalc();
    return this;
  }
  /**
   * Set the current offset to apply to the `points` before applying the `rotation` rotation.
   * @param {Vector} offset The new offset vector
   * @return {Polygon} This for chaining
   */
  setOffset(offset) {
    this.offset = offset;
    this._recalc();
    return this;
  }
  /**Rotates this polygon counter-clockwise around the origin of *its local coordinate system* (i.e. `pos`).
   * Note: This changes the **original** points (so any `rotation` will be applied on top of this rotation).
   * @param {Number} rotation The rotation to rotate (in radians)
   * @return {Polygon} This for chaining
   */
  rotate(rotation) {
    var points = this.points;
    for (var i = 0, len = points.length; i < len; i++) {
      points[i].rotate(rotation);
    }
    this._recalc();
    return this;
  }
  /**
   * Translates the points of this polygon by a specified amount relative to the origin of *its own coordinate
   * system* (i.e. `body.position`)
   * This is most useful to change the "center point" of a polygon. If you just want to move the whole polygon, change
   * the coordinates of `body.position`.
   * Note: This changes the **original** points (so any `offset` will be applied on top of this translation)
   * @param {Number} x The horizontal amount to translate
   * @param {Number} y The vertical amount to translate
   * @return {Polygon} This for chaining
   */
  translate(x, y) {
    var points = this.points;
    for (var i = 0, len = points.length; i < len; i++) {
      points[i].x += x;
      points[i].y += y;
    }
    this._recalc();
    return this;
  }
  /**
   * Computes the calculated collision polygon. Applies the `rotation` and `offset` to the original points then recalculates the
   * edges and normals of the collision polygon.
   * @return {Polygon} This for chaining
   */
  _recalc() {
    // Calculated points - this is what is used for underlying collisions and takes into account
    // the rotation/offset set on the polygon.
    var calcPoints = this.calcPoints;
    // The edges here are the direction of the `n`th edge of the polygon, relative to
    // the `n`th point. If you want to draw a given edge from the edge value, you must
    // first translate to the position of the starting point.
    var edges = this.edges;
    // The normals here are the direction of the normal for the `n`th edge of the polygon, relative
    // to the position of the `n`th point. If you want to draw an edge normal, you must first
    // translate to the position of the starting point.
    var normals = this.normals;
    // Copy the original points array and apply the offset/rotation
    var points = this.points;
    var offset = this.offset;
    var rotation = this._rotation;
    var len = points.length;
    var i;
    for (i = 0; i < len; i++) {
      var calcPoint = calcPoints[i].copy(points[i]);
      calcPoint.x += offset.x;
      calcPoint.y += offset.y;
      if (rotation !== 0) {
        calcPoint.rotate(rotation);
      }
    }
    // Calculate the edges/normals
    for (i = 0; i < len; i++) {
      var p1 = calcPoints[i];
      var p2 = i < len - 1 ? calcPoints[i + 1] : calcPoints[0];
      var e = edges[i].copy(p2).subtract(p1);
      normals[i].copy(e).perp().normalize();
    }
    return this;
  }
};
Object.defineProperty(Polygon.prototype, 'rotation', {
  get: function() {
    return this._rotation;
  },
  set: function(rotation) {
    this._rotation = rotation;
    this._recalc();
  },
});

Object.assign(Rectangle.prototype, {
  rotation: 0,
  /**
   * Returns a polygon whose edges are the same as this rectangle.
   * @return {Polygon} A new Polygon that represents this box
   */
  toPolygon: function toPolygon() {
    var halfWidth = this.width * 0.5;
    var halfHeight = this.height * 0.5;
    return new Polygon([
      new Vector(-halfWidth, -halfHeight), new Vector(halfWidth, -halfHeight),
      new Vector(halfWidth, halfHeight), new Vector(-halfWidth, halfHeight)
    ]);
  },
});

Object.assign(Circle.prototype, {
  rotation: 0,
});

/**
 * Response
 * An object representing the result of an intersection. Contains:
 * - The two objects participating in the intersection
 * - The vector representing the minimum change necessary to extract the first object
 *   from the second one (as well as a unit vector in that direction and the magnitude
 *   of the overlap)
 * - Whether the first object is entirely inside the second, and vice versa.
 */
class Response {
  constructor() {
    this.a = null;
    this.b = null;
    this.aInB = true;
    this.bInA = true;
    this.overlap = Number.MAX_VALUE;
    this.overlapN = new Vector();
    this.overlapV = new Vector();
  }
  /**
   * Set some values of the response back to their defaults.  Call this between tests if
   * you are going to reuse a single Response object for multiple intersection tests (recommented
   * as it will avoid allcating extra memory)
   * @return {Response} This for chaining
   */
  clear() {
    this.aInB = true;
    this.bInA = true;
    this.overlap = Number.MAX_VALUE;
    return this;
  }
};

/**
 * SAT based collision solver
 */
class SATSolver {
  /**
   * Hit test a versus b.
   * @method hitTest
   * @param {Body} a
   * @param {Body} b
   * @return {Boolean} return true, if bodies hit.
   */
  hitTest(a, b, response) {
    // Polygon vs polygon
    if (a.shape.points && b.shape.points) {
      return testPolygonPolygon(a, b, response);
    }
    if (a.shape.radius && b.shape.radius) {
      return testCircleCircle(a, b, response);
    }
    if (a.shape.points && b.shape.radius) {
      return testPolygonCircle(a, b, response);
    }
    else if (a.shape.radius && b.shape.points) {
      return testCirclePolygon(a, b, response);
    }

    throw 'Hit test should not go so far!';
    return false;
  }
  hitResponse(a, b, AvsB, BvsA, response) {
    // Make sure a and b are not reversed
    var uniqueA = (a === response.a ? a : b),
      uniqueB = (b === response.b ? b : a);
    var responseToA = false,
      responseToB = false;
    // Check to see which one or two finally get the response
    if (AvsB && !BvsA) {
      responseToA = uniqueA.collide(uniqueB, response);
    }
    else if (!AvsB && BvsA) {
      responseToB = uniqueB.collide(uniqueA, response);
    }
    else if (AvsB && BvsA) {
      responseToA = uniqueA.collide(uniqueB, response);
      responseToB = uniqueB.collide(uniqueA, response);
    }

    // Only apply response to A if it wants to
    if (responseToA && !responseToB) {
      uniqueA.position.subtract(response.overlapV);
      uniqueA.afterCollide(uniqueB);
    }
    // Only apply response to B if it wants to
    else if (!responseToA && responseToB) {
      uniqueB.position.add(response.overlapV);
      uniqueB.afterCollide(uniqueA);
    }
    // Apply response to both A and B
    else if (responseToA && responseToB) {
      response.overlapV.multiply(0.5);
      uniqueA.position.subtract(response.overlapV);
      uniqueB.position.add(response.overlapV);
      uniqueA.afterCollide(uniqueB);
      uniqueB.afterCollide(uniqueA);
    }
  }
};

class World {
  constructor(x, y) {
    x = typeof x === 'number' ? x : 0;
    y = typeof y === 'number' ? y : 980;
    this.gravity = new Vector(x, y);
    this.bodies = [];
    this.solver = new SATSolver();

    /**
     * SpatialGrid for Broad-Phase Collision.
     * @property {Object} spatialGrid
     */
    this.spatialGrid = new SpatialGrid(this.cellSize, this.bodies);

    /**
     * Cell size of spatial grids
     * @type {Number}
     */
    this.cellSize = 64;
  }

  /**
   * Add body to world.
   * @method addBody
   * @param {Body} body
   */
  addBody(body) {
    body.world = this;
    body._remove = false;
    this.bodies.push(body);
  }

  /**
   * Remove body from world.
   * @method removeBody
   *  @param {Body} body
   */
  removeBody(body) {
    if (!body.world) return;
    body.world = null;
    body._remove = true;
  }

  preUpdate() {}

  /**
   * Update physics world.
   *  @method update
   */
  update(delta) {
    var i, j, bodies = this.spatialGrid.bodies;
    for (i = bodies.length - 1; i >= 0; i--) {
      if (bodies[i]._remove) {
        bodies.splice(i, 1);
      }
      else {
        bodies[i].update(delta);
      }
    }
    this.spatialGrid.update();
    this.spatialGrid.handleCollision(this.solver, this.bodyShouldCollideWithBody);
  }
  bodyShouldCollideWithBody(bodyA, bodyB) {
    return (bodyA.collideAgainst.indexOf(bodyB.collisionGroup) !== -1);
  }
};

class Body extends BaseBody {
  constructor(settings) {
    super(settings);
    this._id = Body.uid++;

    if (this.shape) {
      this.shape.body = this;

      // Convert Rectangle to Polygon
      // TODO: remove this after Rectangle vs Rectanvle completed
      if (this.shape.width) {
        this.shape = this.shape.toPolygon();
      }
    }
  }
  addShape(shape) {
    this.shape = shape;
    // Convert Rectangle to Polygon
    // TODO: remove this after Rectangle vs Rectanvle completed
    if (this.shape.width) {
      this.shape = this.shape.toPolygon();
    }

    shape.body = this;
    return this;
  }
  /**
    Set new collision group for body.
    @method setCollisionGroup
    @param {Number} group
  **/
  setCollisionGroup(group) {
    this.collisionGroup = group;
  }
};
Body.uid = 0;

Object.defineProperty(Body.prototype, 'rotation', {
  get: function() {
    return this.shape ? this.shape.rotation : 0;
  },
  set: function(rotation) {
    this.shape && (this.shape.rotation = rotation);
  },
});

class SpatialGrid {
  constructor(cellSize, bodies) {
    /**
     * Grid cell size in pixel
     * @type {Number}
     */
    this.pxCellSize = cellSize || 64;
    /**
     * Body list
     * @type {Array}
     */
    this.bodies = bodies;
    /**
     * Cell matrix
     * @type {Array<Array>}
     */
    this.grid = [[]];

    /**
     * Response object for testing reuse
     * @type {Response}
     */
    this.response = new Response();

    /**
     * How many collision tests occured within current step
     * @type {Number}
     */
    this.collisionTests = 0;
    /**
     * How many cells allocated for collision test
     * @type {Number}
     */
    this.allocatedCells = 0;
    /**
     * How many collision checking occured within current step
     * @type {Number}
     */
    this.hashChecks = 0;
    this.totalCells = 0;
  }

  update() {
    var cXEntityMin, cXEntityMax, cYEntityMin, cYEntityMax, i, j, body, cX, cY, gridCol, gridCell;

    this.allocatedCells = 0;

    // Construct grid
    // NOTE: this is a purposeful use of the Array() constructor
    this.grid = Object.create(null);

    // Insert all bodies into grid
    for (i = 0; i < this.bodies.length; i++) {
      body = this.bodies[i];

      // Find extremes of cells that body overlaps
      // Subtract min to shift grid to avoid negative numbers
      var top, right, bottom, left;
      if (body.shape.points) {
        var axis = T_VECTORS.pop(),
          res = T_ARRAYS.pop();
        flattenPointsOn(body.shape.points, axis.set(1, 0), res);
        left = res[0] + body.position.x;
        right = res[1] + body.position.x;
        flattenPointsOn(body.shape.points, axis.set(0, 1), res);
        top = res[0] + body.position.y;
        bottom = res[1] + body.position.y;
        T_VECTORS.push(axis);
        T_ARRAYS.push(res);
      }
      else if (body.shape.radius) {
        left = body.position.x - body.shape.radius;
        right = body.position.x + body.shape.radius;
        top = body.position.y - body.shape.radius;
        bottom = body.position.y + body.shape.radius;
      }
      cXEntityMin = Math.floor(left / this.pxCellSize);
      cXEntityMax = Math.floor(right / this.pxCellSize);
      cYEntityMin = Math.floor(top / this.pxCellSize);
      cYEntityMax = Math.floor(bottom / this.pxCellSize);

      // Insert body into each cell it overlaps
      // We're looping to make sure that all cells between extremes are found
      for (cX = cXEntityMin; cX <= cXEntityMax; cX++) {
        // Make sure a column exists, initialize if not
        gridCol = this.grid[cX];
        if (!gridCol) {
          this.grid[cX] = gridCol = Object.create(null);
        }

        // Loop through each cell in this column
        for (cY = cYEntityMin; cY <= cYEntityMax; cY++) {
          // Ensure we have a bucket to put bodies into for this cell
          gridCell = gridCol[cY];
          if (!gridCell) {
            gridCol[cY] = gridCell = [];

            // This is for stats purposes only
            this.allocatedCells += 1;
          }

          // Add body to cell
          gridCell.push(body);
        }
      }
    }
  }

  handleCollision(solver, shouldCollide) {
    var checked = Object.create(null),
      bodyA, bodyB, hashA, hashB, i, j, k, l, coIdx, gridCol, gridCell,
      aShouldCollideWithB = false, bShouldCollideWithA = false;

    // Reset counts, for debug/comparison purposes
    this.collisionTests = 0;
    this.hashChecks = 0;

    // For every column in the grid...
    for (i in this.grid) {
      // Note: these lines are disabled since grid is now a simple "hasn"
      // if (!this.grid.hasOwnProperty(i)) {
      //     continue;
      // }

      gridCol = this.grid[i];

      // For every cell within a column of the grid...
      for (j in gridCol) {
        // Note: these lines are disabled since gridCol is now a simple "hasn"
        // if (!gridCol.hasOwnProperty(j)) {
        //     continue;
        // }
        gridCell = gridCol[j];

        // For every object in a cell...
        for (k = 0; k < gridCell.length; k++) {
          bodyA = gridCell[k];
          if (bodyA._remove) continue;

          // For every other object in a cell...
          for (l = k + 1; l < gridCell.length; l++) {
            bodyB = gridCell[l];
            if (bodyB._remove) continue;

            // Skip if they should not collide with each other
            aShouldCollideWithB = shouldCollide(bodyA, bodyB);
            bShouldCollideWithA = shouldCollide(bodyB, bodyA);
            if (!aShouldCollideWithB && !bShouldCollideWithA) {
              continue;
            }

            // Create a unique key to mark this pair
            hashA = bodyA._id + ':' + bodyB._id;
            hashB = bodyB._id + ':' + bodyA._id;

            this.hashChecks += 2;

            if (!checked[hashA] && !checked[hashB]) {
              checked[hashA] = checked[hashB] = true;

              this.collisionTests += 1;

              // Test collision
              this.response.clear();
              if (solver.hitTest(bodyA, bodyB, this.response)) {
                // Solve collision response
                solver.hitResponse(bodyA, bodyB, aShouldCollideWithB, bShouldCollideWithA, this.response);
              }
            }
          }
        }
      }
    }
  }
};

// Helper Functions ------------------------------------

/**
 * Flattens the specified array of points onto a unit vector axis,
 * resulting in a one dimensional range of the minimum and
 * maximum value on that axis.
 * @param {Array<Vector>} points The points to flatten
 * @param {Vector} normal The unit vector axis to flatten on
 * @param {Array<Number>} result An array.  After calling this function,
 *   result[0] will be the minimum value,
 *   result[1] will be the maximum value
 */
function flattenPointsOn(points, normal, result) {
  var min = Number.MAX_VALUE;
  var max = -Number.MAX_VALUE;
  var len = points.length;
  for (var i = 0; i < len; i++) {
    // The magnitude of the projection of the point onto the normal
    var dot = points[i].dot(normal);
    if (dot < min) { min = dot; }
    if (dot > max) { max = dot; }
  }
  result[0] = min; result[1] = max;
}

/**
 * Check whether two convex polygons are separated by the specified
 * axis (must be a unit vector).
 * @param {Vector} aPos The position of the first polygon
 * @param {Vector} bPos The position of the second polygon
 * @param {Array<Vector>} aPoints The points in the first polygon
 * @param {Array<Vector>} bPoints The points in the second polygon
 * @param {Vector} axis The axis (unit sized) to test against. The points of both polygons
 *   will be projected onto this axis
 * @param {Response=} response A Response object (optional) which will be populated
 *   if the axis is not a separating axis
 * @return {Boolean} true if it is a separating axis, false otherwise.  If false,
 *   and a response is passed in, information about how much overlap and
 *   the direction of the overlap will be populated
 */
function isSeparatingAxis(aPos, bPos, aPoints, bPoints, axis, response) {
  var rangeA = T_ARRAYS.pop();
  var rangeB = T_ARRAYS.pop();
  // The magnitude of the offset between the two polygons
  var offsetV = T_VECTORS.pop().copy(bPos).subtract(aPos);
  var projectedOffset = offsetV.dot(axis);
  // Project the polygons onto the axis.
  flattenPointsOn(aPoints, axis, rangeA);
  flattenPointsOn(bPoints, axis, rangeB);
  // Move B's range to its position relative to A.
  rangeB[0] += projectedOffset;
  rangeB[1] += projectedOffset;
  // Check if there is a gap. If there is, this is a separating axis and we can stop
  if (rangeA[0] > rangeB[1] || rangeB[0] > rangeA[1]) {
    T_VECTORS.push(offsetV);
    T_ARRAYS.push(rangeA);
    T_ARRAYS.push(rangeB);
    return true;
  }
  // This is not a separating axis. If we're calculating a response, calculate the overlap.
  if (response) {
    var overlap = 0;
    // A starts further left than B
    if (rangeA[0] < rangeB[0]) {
      response.aInB = false;
      // A ends before B does. We have to pull A out of B
      if (rangeA[1] < rangeB[1]) {
        overlap = rangeA[1] - rangeB[0];
        response.bInA = false;
        // B is fully inside A.  Pick the shortest way out.
      }
      else {
        var option1 = rangeA[1] - rangeB[0];
        var option2 = rangeB[1] - rangeA[0];
        overlap = option1 < option2 ? option1 : -option2;
      }
    // B starts further left than A
    }
    else {
      response.bInA = false;
      // B ends before A ends. We have to push A out of B
      if (rangeA[1] > rangeB[1]) {
        overlap = rangeA[0] - rangeB[1];
        response.aInB = false;
      // A is fully inside B.  Pick the shortest way out.
      }
      else {
        var option1 = rangeA[1] - rangeB[0];
        var option2 = rangeB[1] - rangeA[0];
        overlap = option1 < option2 ? option1 : -option2;
      }
    }
    // If this is the smallest amount of overlap we've seen so far, set it as the minimum overlap.
    var absOverlap = Math.abs(overlap);
    if (absOverlap < response.overlap) {
      response.overlap = absOverlap;
      response.overlapN.copy(axis);
      if (overlap < 0) {
        response.overlapN.reverse();
      }
    }
  }
  T_VECTORS.push(offsetV);
  T_ARRAYS.push(rangeA);
  T_ARRAYS.push(rangeB);
  return false;
}

/**
 * Calculates which Vornoi region a point is on a line segment.
 * It is assumed that both the line and the point are relative to `(0,0)`
 *            |       (0)      |
 *     (-1)  [S]--------------[E]  (1)
 *            |       (0)      |
 * @param {Vector} line The line segment
 * @param {Vector} point The point
 * @return {number} LEFT_VORNOI_REGION (-1) if it is the left region,
 *         MIDDLE_VORNOI_REGION (0) if it is the middle region,
 *         RIGHT_VORNOI_REGION (1) if it is the right region
 */
function vornoiRegion(line, point) {
  var len2 = line.squaredLength();
  var dp = point.dot(line);
  // If the point is beyond the start of the line, it is in the
  // left vornoi region.
  if (dp < 0) { return LEFT_VORNOI_REGION; }
  // If the point is beyond the end of the line, it is in the
  // right vornoi region.
  else if (dp > len2) { return RIGHT_VORNOI_REGION; }
  // Otherwise, it's in the middle one.
  else { return MIDDLE_VORNOI_REGION; }
}
// Constants for Vornoi regions
var LEFT_VORNOI_REGION = -1;
var MIDDLE_VORNOI_REGION = 0;
var RIGHT_VORNOI_REGION = 1;

// Collision Tests ---------------------------------------

/**
 * Check if a point is inside a circle.
 * @param {Vector} p The point to test
 * @param {Circle} c The circle to test
 * @return {Boolean} true if the point is inside the circle, false if it is not
 */
function pointInCircle(p, c) {
  var differenceV = T_VECTORS.pop().copy(p).subtract(c.position);
  var radiusSq = c.shape.radius * c.shape.radius;
  var distanceSq = differenceV.squaredLength();
  T_VECTORS.push(differenceV);
  // If the distance between is smaller than the radius then the point is inside the circle.
  return distanceSq <= radiusSq;
}

/**
 * Check if a point is inside a convex polygon.
 * @param {Vector} p The point to test
 * @param {Polygon} poly The polygon to test
 * @return {Boolean} true if the point is inside the polygon, false if it is not
 */
function pointInPolygon(p, poly) {
  UNIT_SQUARE.position.copy(p);
  T_RESPONSE.clear();
  var result = testPolygonPolygon(UNIT_SQUARE, poly, T_RESPONSE);
  if (result) {
    result = T_RESPONSE.aInB;
  }
  return result;
}

/**
 * Check if two circles collide.
 * @param {Body} a The first circle body
 * @param {Body} b The second circle body
 * @param {Response=} response Response object (optional) that will be populated if
 *   the circles intersect
 * @return {Boolean} true if the circles intersect, false if they don't
 */
function testCircleCircle(a, b, response) {
  // Check if the distance between the centers of the two
  // circles is greater than their combined radius.
  var differenceV = T_VECTORS.pop().copy(b.position).subtract(a.position);
  var totalRadius = a.shape.radius + b.shape.radius;
  var totalRadiusSq = totalRadius * totalRadius;
  var distanceSq = differenceV.squaredLength();
  // If the distance is bigger than the combined radius, they don't intersect.
  if (distanceSq > totalRadiusSq) {
    T_VECTORS.push(differenceV);
    return false;
  }
  // They intersect.  If we're calculating a response, calculate the overlap.
  if (response) {
    var dist = Math.sqrt(distanceSq);
    response.a = a;
    response.b = b;
    response.overlap = totalRadius - dist;
    response.overlapN.copy(differenceV.normalize());
    response.overlapV.copy(differenceV).multiply(response.overlap);
    response.aInB = a.shape.radius <= b.shape.radius && dist <= b.shape.radius - a.shape.radius;
    response.bInA = b.shape.radius <= a.shape.radius && dist <= a.shape.radius - b.shape.radius;
  }
  T_VECTORS.push(differenceV);
  return true;
}

/**
 * Check if a polygon and a circle collide.
 * @param {Polygon} polygon The polygon
 * @param {Circle} circle The circle
 * @param {Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
function testPolygonCircle(polygon, circle, response) {
  // Get the position of the circle relative to the polygon.
  var circlePos = T_VECTORS.pop().copy(circle.position).subtract(polygon.position);
  var radius = circle.shape.radius;
  var radius2 = radius * radius;
  var points = polygon.shape.calcPoints;
  var len = points.length;
  var edge = T_VECTORS.pop();
  var point = T_VECTORS.pop();

  // For each edge in the polygon:
  for (var i = 0; i < len; i++) {
    var next = i === len - 1 ? 0 : i + 1;
    var prev = i === 0 ? len - 1 : i - 1;
    var overlap = 0;
    var overlapN = null;

    // Get the edge.
    edge.copy(polygon.shape.edges[i]);
    // Calculate the center of the circle relative to the starting point of the edge.
    point.copy(circlePos).subtract(points[i]);

    // If the distance between the center of the circle and the point
    // is bigger than the radius, the polygon is definitely not fully in
    // the circle.
    if (response && point.squaredLength() > radius2) {
      response.aInB = false;
    }

    // Calculate which Vornoi region the center of the circle is in.
    var region = vornoiRegion(edge, point);
    // If it's the left region:
    if (region === LEFT_VORNOI_REGION) {
      // We need to make sure we're in the RIGHT_VORNOI_REGION of the previous edge.
      edge.copy(polygon.shape.edges[prev]);
      // Calculate the center of the circle relative the starting point of the previous edge
      var point2 = T_VECTORS.pop().copy(circlePos).subtract(points[prev]);
      region = vornoiRegion(edge, point2);
      if (region === RIGHT_VORNOI_REGION) {
        // It's in the region we want.  Check if the circle intersects the point.
        var dist = point.length();
        if (dist > radius) {
          // No intersection
          T_VECTORS.push(circlePos);
          T_VECTORS.push(edge);
          T_VECTORS.push(point);
          T_VECTORS.push(point2);
          return false;
        }
        else if (response) {
          // It intersects, calculate the overlap.
          response.bInA = false;
          overlapN = point.normalize();
          overlap = radius - dist;
        }
      }
      T_VECTORS.push(point2);
    }
    // If it's the right region:
    else if (region === RIGHT_VORNOI_REGION) {
      // We need to make sure we're in the left region on the next edge
      edge.copy(polygon.shape.edges[next]);
      // Calculate the center of the circle relative to the starting point of the next edge.
      point.copy(circlePos).subtract(points[next]);
      region = vornoiRegion(edge, point);
      if (region === LEFT_VORNOI_REGION) {
        // It's in the region we want.  Check if the circle intersects the point.
        var dist = point.length();
        if (dist > radius) {
          // No intersection
          T_VECTORS.push(circlePos);
          T_VECTORS.push(edge);
          T_VECTORS.push(point);
          return false;
        }
        else if (response) {
          // It intersects, calculate the overlap.
          response.bInA = false;
          overlapN = point.normalize();
          overlap = radius - dist;
        }
      }
    }
    // Otherwise, it's the middle region:
    else {
      // Need to check if the circle is intersecting the edge,
      // Change the edge into its "edge normal".
      var normal = edge.perp().normalize();
      // Find the perpendicular distance between the center of the
      // circle and the edge.
      var dist = point.dot(normal);
      var distAbs = Math.abs(dist);
      // If the circle is on the outside of the edge, there is no intersection.
      if (dist > 0 && distAbs > radius) {
        // No intersection
        T_VECTORS.push(circlePos);
        T_VECTORS.push(normal);
        T_VECTORS.push(point);
        return false;
      }
      else if (response) {
        // It intersects, calculate the overlap.
        overlapN = normal;
        overlap = radius - dist;
        // If the center of the circle is on the outside of the edge, or part of the
        // circle is on the outside, the circle is not fully inside the polygon.
        if (dist >= 0 || overlap < 2 * radius) {
          response.bInA = false;
        }
      }
    }

    // If this is the smallest overlap we've seen, keep it.
    // (overlapN may be null if the circle was in the wrong Vornoi region).
    if (overlapN && response && Math.abs(overlap) < Math.abs(response.overlap)) {
      response.overlap = overlap;
      response.overlapN.copy(overlapN);
    }
  }

  // Calculate the final overlap vector - based on the smallest overlap.
  if (response) {
    response.a = polygon;
    response.b = circle;
    response.overlapV.copy(response.overlapN).multiply(response.overlap);
  }
  T_VECTORS.push(circlePos);
  T_VECTORS.push(edge);
  T_VECTORS.push(point);
  return true;
}

/**
 * Check if a circle and a polygon collide.
 *
 * **NOTE:** This is slightly less efficient than polygonCircle as it just
 * runs polygonCircle and reverses everything at the end.
 *
 * @param {Circle} circle The circle
 * @param {Polygon} polygon The polygon
 * @param {Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
function testCirclePolygon(circle, polygon, response) {
  // Test the polygon against the circle.
  var result = testPolygonCircle(polygon, circle, response);
  if (result && response) {
    // Swap A and B in the response.
    var a = response.a;
    var aInB = response.aInB;
    response.overlapN.reverse();
    response.overlapV.reverse();
    response.a = response.b;
    response.b = a;
    response.aInB = response.bInA;
    response.bInA = aInB;
  }
  return result;
}

/**
 * Checks whether polygons collide.
 * @param {Polygon} a The first polygon
 * @param {Polygon} b The second polygon
 * @param {Response=} response Response object (optional) that will be populated if
 *   they interset
 * @return {Boolean} true if they intersect, false if they don't
 */
function testPolygonPolygon(a, b, response) {
  var aPoints = a.shape.calcPoints;
  var aLen = aPoints.length;
  var bPoints = b.shape.calcPoints;
  var bLen = bPoints.length;
  // If any of the edge normals of A is a separating axis, no intersection.
  for (var i = 0; i < aLen; i++) {
    if (isSeparatingAxis(a.position, b.position, aPoints, bPoints, a.shape.normals[i], response)) {
      return false;
    }
  }
  // If any of the edge normals of B is a separating axis, no intersection.
  for (var i = 0;i < bLen; i++) {
    if (isSeparatingAxis(a.position, b.position, aPoints, bPoints, b.shape.normals[i], response)) {
      return false;
    }
  }
  // Since none of the edge normals of A or B are a separating axis, there is an intersection
  // and we've already calculated the smallest overlap (in isSeparatingAxis).  Calculate the
  // final overlap vector.
  if (response) {
    response.a = a;
    response.b = b;
    response.overlapV.copy(response.overlapN).multiply(response.overlap);
  }
  return true;
}

// Object Pools -----------------------------------------

/**
 * A pool of `Vector` objects that are used in calculations to avoid
 * allocating memory.
 * @type {Array<Vector>}
 */
var T_VECTORS = [];
for (var i = 0; i < 10; i++) { T_VECTORS.push(new Vector()); }

/**
 * A pool of arrays of numbers used in calculations to avoid allocating
 * memory.
 * @type {Array<Array<Bumber>>}
 */
var T_ARRAYS = [];
for (var i = 0; i < 5; i++) { T_ARRAYS.push([]); }

/**
 * Temporary response used for polygon hit detection.
 * @type {Response}
 */
var T_RESPONSE = new Response();

/**
 * Unit square polygon used for polygon hit detection.
 * @type {Polygon}
 */
var UNIT_SQUARE = new Body({
  position: new Vector(),
  shape: new Rectangle(1, 1).toPolygon(),
});

// Register as sub-system
Scene.registerSystem('Physics', {
  init: function init(scene) {
    scene.world = new World();
  },
  preUpdate: function preUpdate(scene, delta) {
    scene.world.preUpdate(delta * 0.001);
  },
  update: function update(scene, delta) {
    scene.world.update(delta * 0.001);
  },
});

export default {
  // Shapes
  Rectangle,
  Circle,
  Polygon,

  // Core objects
  Body,
  World,

  // Core of collision detection
  SATSolver,
  SpatialGrid,
  Response,

  // Overlap testing functions
  testCircleCircle,
  testPolygonCircle,
  testCirclePolygon,
  testPolygonPolygon,

  // Utils
  flattenPointsOn,
  isSeparatingAxis,
  vornoiRegion,
  pointInCircle,
  pointInPolygon,
};
