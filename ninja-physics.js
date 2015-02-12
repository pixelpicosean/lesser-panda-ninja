game.module(
    'plugins.ninja-physics'
)
.require(
    'engine.physics'
)
.body(function() { 'use strict';

    game.Vector.inject({
        /**
         * Change this vector to be perpendicular to what it was before. (Effectively
         * roatates it 90 degrees in a clockwise direction)
         * @return {game.Vector} This for chaining
         */
        perp: function() {
            var x = this.x;
            this.x = this.y;
            this.y = -x;
            return this;
        },
        /**
         * Reverse this vector
         * @return {game.Vector} This for chaining
         */
        reverse: function() {
            this.x = -this.x;
            this.y = -this.y;
            return this;
        },
        /**
         * Scale this vector. An independant scaling factor can be provided
         * for each axis, or a single scaling factor that will scale both `x` and `y`.
         * @param  {Number} x The scaling factor in the x direction
         * @param  {Number} y The scaling factor in the y direction
         * @return {game.Vector} This for chaining
         */
        scale: function(x, y) {
            this.x *= x;
            this.y *= y || x;
            return this;
        },
        /**
         * Project this vector on to another vector.
         * @param {game.Vector} other The vector to project onto
         * @return {game.Vector} This for chaining
         */
        project: function(other) {
            var amt = this.dot(other) / other.squaredLength();
            this.x = amt * other.x;
            this.y = amt * other.y;
            return this;
        },
        /**
         * Project this vector onto a vector of unit length. This is slightly more efficient
         * than `project` when dealing with unit vectors.
         * @param {game.Vector} other The unit vector to project onto
         * @return {game.Vector} This for chaining
         */
        projectN: function(other) {
            var amt = this.dot(other);
            this.x = amt * other.x;
            this.y = amt * other.y;
            return this;
        },
        /**
         * Reflect this vector on an arbitrary axis.
         * @param {game.Vector} axis The vector representing the axis
         * @return {game.Vector} This for chaining
         */
        reflect: function(axis) {
            var x = this.x,
                y = this.y;
            this.project(axis).scale(2);
            this.x -= x;
            this.y -= y;
            return this;
        },
        /**
         * Reflect this vector on an arbitrary axis (represented by a unit vector). This is
         * slightly more efficient than `reflect` when dealing with an axis that is a unit vector.
         * @param {game.Vector} axis The unit vector representing the axis
         * @return {game.Vector} This for chaining
         */
        reflectN: function(axis) {
            var x = this.x,
                y = this.y;
            this.projectN(axis).scale(2);
            this.x -= x;
            this.y -= y;
            return this;
        },
        /**
         * Get the squared length of this vector.
         * @return {Number} The length^2 of this vector
         */
        squaredLength: function() {
            return this.dot();
        }
    });

    // TODO: finish polygon related methods
    game.createClass('Polygon', {
        points: [],
        calcPoints: [],
        edges: [],
        normals: [],
        offset: null,
        angle: 0,
        init: function(points) {
            this.offset = new game.Vector();
            this.setPoints(points || []);
        },
        /**
         * Set the points of the polygon.
         * @param {Array<game.Vector>=} points An array of vectors representing the points in the polygon,
         *   in counter-clockwise order
         * @return {game.Polygon} This for chaining
         */
        setPoints: function(points) {
            // Only re-allocate if this is a new polygon or the number of points has changed.
            var lengthChanged = !this.points || this.points.length !== points.length;
            if (lengthChanged) {
                var calcPoints = this.calcPoints = [];
                var edges = this.edges = [];
                var normals = this.normals = [];
                // Allocate the vector arrays for the calculated properties
                for (var i = 0, len = points.length; i < len; i++) {
                    calcPoints.push(new game.Vector());
                    edges.push(new game.Vector());
                    normals.push(new game.Vector());
                }
            }
            this.points = points;
            this._recalc();
            return this;
        },
        /**
         * Set the current rotation angle of the polygon.
         * @param {Number} angle The current rotation angle (in radians)
         * @return {game.Polygon} This for chaining
         */
        setAngle: function(angle) {
            this.angle = angle;
            this._recalc();
            return this;
        },
        /**
         * Set the current offset to apply to the `points` before applying the `angle` rotation.
         * @param {game.Vector} offset The new offset vector
         * @return {game.Polygon} This for chaining
         */
        setOffset: function(offset) {
            this.offset = offset;
            this._recalc();
            return this;
        },
        /**Rotates this polygon counter-clockwise around the origin of *its local coordinate system* (i.e. `pos`).
         * Note: This changes the **original** points (so any `angle` will be applied on top of this rotation).
         * @param {Number} angle The angle to rotate (in radians)
         * @return {game.Polygon} This for chaining
         */
        rotate: function(angle) {
            var points = this.points;
            for (var i = 0, len = points.length; i < len; i++) {
                points[i].rotate(angle);
            }
            this._recalc();
            return this;
        },
        /**
         * Translates the points of this polygon by a specified amount relative to the origin of *its own coordinate
         * system* (i.e. `body.position`)
         * This is most useful to change the "center point" of a polygon. If you just want to move the whole polygon, change
         * the coordinates of `body.position`.
         * Note: This changes the **original** points (so any `offset` will be applied on top of this translation)
         * @param {Number} x The horizontal amount to translate
         * @param {Number} y The vertical amount to translate
         * @return {game.Polygon} This for chaining
         */
        translate: function(x, y) {
            var points = this.points;
            for (var i = 0, len = points.length; i < len; i++) {
                points[i].x += x;
                points[i].y += y;
            }
            this._recalc();
            return this;
        },
        /**
         * Computes the calculated collision polygon. Applies the `angle` and `offset` to the original points then recalculates the
         * edges and normals of the collision polygon.
         * @return {game.Polygon} This for chaining
         */
        _recalc: function() {
            // Calculated points - this is what is used for underlying collisions and takes into account
            // the angle/offset set on the polygon.
            var calcPoints = this.calcPoints;
            // The edges here are the direction of the `n`th edge of the polygon, relative to
            // the `n`th point. If you want to draw a given edge from the edge value, you must
            // first translate to the position of the starting point.
            var edges = this.edges;
            // The normals here are the direction of the normal for the `n`th edge of the polygon, relative
            // to the position of the `n`th point. If you want to draw an edge normal, you must first
            // translate to the position of the starting point.
            var normals = this.normals;
            // Copy the original points array and apply the offset/angle
            var points = this.points;
            var offset = this.offset;
            var angle = this.angle;
            var len = points.length;
            var i;
            for (i = 0; i < len; i++) {
                var calcPoint = calcPoints[i].copy(points[i]);
                calcPoint.x += offset.x;
                calcPoint.y += offset.y;
                if (angle !== 0) {
                    calcPoint.rotate(angle);
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
    });

    game.Rectangle.inject({
        /**
         * Returns a polygon whose edges are the same as this rectangle.
         * @return {game.Polygon} A new Polygon that represents this box
         */
        toPolygon: function() {
            return new game.Polygon([
                new game.Vector(), new game.Vector(this.width, 0),
                new game.Vector(this.width, this.height), new game.Vector(0, this.height)
            ]);
        }
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
    game.createClass('Response', {
        a: null,
        b: null,
        overlapN: null,
        overlapV: null,
        aInB: true,
        bInA: true,
        overlap: Number.MAX_VALUE,
        init: function() {
            this.overlapN = new game.Vector();
            this.overlapV = new game.Vector();
        },
        /**
         * Set some values of the response back to their defaults.  Call this between tests if
         * you are going to reuse a single Response object for multiple intersection tests (recommented
         * as it will avoid allcating extra memory)
         * @return {game.Response} This for chaining
         */
        clear: function() {
            this.aInB = true;
            this.bInA = true;
            this.overlap = Number.MAX_VALUE;
            return this;
        }
    });

    game.World.inject({
        /**
         * SpatialGrid for Broad-Phase Collision.
         * @property {Object} spatialGrid
         */
        spatialGrid: null,
        /**
         * Cell size of spatial grids
         * @type {Number}
         */
        cellSize: 64,

        init: function(x, y) {
            this.gravity = new game.Vector(x || 0, y || 980);
            this.solver = new game.CollisionSolver();
            // Initial size of the grid is 1x1 in cell
            this.spatialGrid = new game.SpatialGrid(this.cellSize);
        },

        /**
         * Add body to world.
         * @method addBody
         * @param {game.Body} body
         */
        addBody: function(body) {
            body.world = this;
            body._remove = false;
            this.spatialGrid.addBody(body);
            if (game.debugDraw && body.shape) {
                game.debugDraw.addBody(body);
            }
        },

        /**
         * Remove body from world.
         * @method removeBody
         *  @param {game.Body} body
         */
        removeBody: function(body) {
            if (!body.world) return;
            body.world = null;
            this.spatialGrid.removeBody(body);
        },

        /**
         * Update physics world.
         *  @method update
         */
        update: function() {
            var i, j,
                bodies = this.spatialGrid.bodies;
            for (i = bodies.length - 1; i >= 0; i--) {
                if (bodies[i]._remove) {
                    bodies.splice(i, 1);
                }
                else {
                    bodies[i].update();
                }
            }
            this.spatialGrid.update();
            this.spatialGrid.handleCollision(this.solver, this.bodyShouldCollideWithBody);
        },
        bodyShouldCollideWithBody: function(bodyA, bodyB) {
            return (bodyA.collideAgainst.indexOf(bodyB.collisionGroup) !== -1);
        }
    });

    game.Body.inject({
        _id: 0,
        init: function(settings) {
            this._super(settings);
            this._id = game.Body.uid++;

            if (this.shape) {
                this.shape.body = this;
            }
        },
        addShape: function(shape) {
            this.shape = shape;
            shape.body = this;
            return this;
        }
    });
    game.Body.uid = 0;

    game.createClass('SpatialGrid', {
        /**
         * Position of the grid left-top point
         * @type {game.Vector}
         */
        min: null,
        /**
         * Position of the grid bottom-right point
         * @type {game.Vector}
         */
        max: null,
        /**
         * Grid cell size in pixel
         * @type {Number}
         */
        pxCellSize: 64,
        /**
         * Body list
         * @type {Array}
         */
        bodies: null,
        /**
         * Cell matrix
         * @type {Array<Array>}
         */
        grid: null,

        /**
         * How many collision tests occured within current step
         * @type {Number}
         */
        collisionTests: 0,
        /**
         * How many cells should be in total for current grid setting
         * @type {Number}
         */
        totalCells: 0,
        /**
         * How many cells allocated for collision test
         * @type {Number}
         */
        allocatedCells: 0,
        /**
         * How many collision checking occured within current step
         * @type {Number}
         */
        hashChecks: 0,

        init: function(cellSize) {
            this.bodies = [];

            this.min = new game.Vector(0, 0);
            this.max = new game.Vector(cellSize, cellSize);
            this.pxCellSize = cellSize;
            this.grid = [[]];

            // These are purely for reporting purposes
            this.collisionTests = 0;
            this.totalCells = 0;
            this.allocatedCells = 0;
            this.hashChecks = 0;
        },
        update: function() {
            var cGridWidth = Math.floor((this.max.x - this.min.x) / this.pxCellSize),
                cGridHeight = Math.floor((this.max.y - this.min.y) / this.pxCellSize),
                cXEntityMin, cXEntityMax, cYEntityMin, cYEntityMax, i, j, body, cX, cY, gridCol, gridCell;

            // The total number of cells this grid will contain.
            this.totalCells = cGridWidth * cGridHeight;
            this.allocatedCells = 0;

            // Construct grid
            // NOTE: this is a purposeful use of the Array() constructor
            this.grid = Array(cGridWidth);

            // Insert all bodies into grid
            for (i = 0; i < this.bodies.length; i++) {
                body = this.bodies[i];

                // If body is outside the grid extents, then ignore it
                if (body.position.x < this.min.x ||
                    body.position.x > this.max.x ||
                    body.position.y < this.min.y ||
                    body.position.y > this.max.y) {
                    continue;
                }

                // Find extremes of cells that body overlaps
                // Subtract min to shift grid to avoid negative numbers
                var bodyWidth, bodyHeight;
                if (body.shape.width) {
                    bodyWidth = body.shape.width;
                    bodyHeight = body.shape.height;
                }
                else if (body.shape.radius) {
                    bodyWidth = bodyHeight = body.shape.radius;
                }
                cXEntityMin = Math.floor((body.position.x - this.min.x) / this.pxCellSize);
                cXEntityMax = Math.floor((body.position.x + bodyWidth - this.min.x) / this.pxCellSize);
                cYEntityMin = Math.floor((body.position.y - this.min.y) / this.pxCellSize);
                cYEntityMax = Math.floor((body.position.y + bodyHeight - this.min.y) / this.pxCellSize);

                // Insert body into each cell it overlaps
                // We're looping to make sure that all cells between extremes are found
                for (cX = cXEntityMin; cX <= cXEntityMax; cX++) {
                    // Make sure a column exists, initialize if not to grid height length
                    // NOTE: again, a purposeful use of the Array constructor
                    if (!this.grid[cX]) {
                        this.grid[cX] = Array(cGridHeight);
                    }

                    gridCol = this.grid[cX];

                    // Loop through each cell in this column
                    for (cY = cYEntityMin; cY <= cYEntityMax; cY++) {
                        // Ensure we have a bucket to put bodies into for this cell
                        if (!gridCol[cY]) {
                            gridCol[cY] = [];

                            // This is for stats purposes only
                            this.allocatedCells += 1;
                        }

                        gridCell = gridCol[cY];

                        // Add body to cell
                        gridCell.push(body);
                    }
                }
            }
        },

        addBody: function(body) {
            var halfWidth = 1, halfHeight = 1;
            if (body.shape.radius) {
                halfWidth = halfHeight = body.shape.radius;
            }
            else {
                if (body.shape.width) {
                    halfWidth = body.shape.width * 0.5;
                }
                if (body.shape.height) {
                    halfHeight = body.shape.height * 0.5;
                }
            }

            // Bounds should be re-initialized when the bodies list is empty
            if (this.bodies.length === 0) {
                this.min.x = Math.floor(body.position.x - halfWidth);
                this.min.y = Math.floor(body.position.y - halfHeight);
                this.max.x = Math.floor(body.position.x + halfWidth);
                this.max.y = Math.floor(body.position.y + halfHeight);
            }
            else {
                // Fast update bounds instead of iterator whole list
                if (body.position.x - halfWidth < this.min.x) {
                    this.min.x = Math.floor(body.position.x - halfWidth);
                }
                else if (body.position.x + halfWidth > this.max.x) {
                    this.max.x = Math.floor(body.position.x + halfWidth);
                }
                if (body.position.y - halfHeight < this.min.y) {
                    this.min.y = Math.floor(body.position.y - halfHeight);
                }
                else if (body.position.y + halfHeight > this.max.y) {
                    this.max.y = Math.floor(body.position.y + halfHeight);
                }
            }

            this.bodies.push(body);
        },

        removeBody: function(body) {
            body._remove = true;
            this.updateBounds();
        },

        handleCollision: function(solver, shouldCollide) {
            var checked = {},
                bodyA, bodyB, hashA, hashB, i, j, k, l, coIdx, gridCol, gridCell,
                aShouldCollideWithB = false, bShouldCollideWithA = false;

            // Reset counts, for debug/comparison purposes
            this.collisionTests = 0;
            this.hashChecks = 0;

            // For every column in the grid...
            for (i = 0; i < this.grid.length; i++) {
                gridCol = this.grid[i];

                // Ignore columns that have no cells
                if (!gridCol) {
                    continue;
                }

                // For every cell within a column of the grid...
                for (j = 0; j < gridCol.length; j++) {
                    gridCell = gridCol[j];

                    // Ignore cells that have no objects
                    if (!gridCell) {
                        continue;
                    }

                    // For every object in a cell...
                    for (k = 0; k < gridCell.length; k++) {
                        bodyA = gridCell[k];

                        // For every other object in a cell...
                        for (l = 0; l < gridCell.length; l++) {
                            bodyB = gridCell[l];

                            if ((bodyA === bodyB) || !shouldCollide(bodyA, bodyB)) {
                                continue;
                            }

                            // Create a unique key to mark this check
                            hashA = bodyA._id + ':' + bodyB._id;

                            this.hashChecks += 1;

                            if (!checked[hashA]) {
                                checked[hashA] = true;

                                this.collisionTests += 1;

                                // Test collision
                                if (solver.hitTest(bodyA, bodyB)) {
                                    // Solve collision
                                    if (solver.hitResponse(bodyA, bodyB)) {
                                        bodyA.afterCollide(bodyB);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        },

        updateBounds: function() {
            if (this.bodies.length === 0) {
                return;
            }

            var i, len, body, halfWidth, halfHeight;
            for (i = 0, len = this.bodies.length; i < len; i++) {
                body = this.bodies[i];
                if (body.shape.width) {
                    halfWidth = body.shape.width * 0.5;
                    halfHeight = body.shape.height * 0.5;
                }
                else if (body.shape.radius) {
                    halfWidth = halfHeight = body.shape.radius;
                }
                this.min.x = Math.min(this.min.x, body.position.x - halfWidth);
                this.min.y = Math.min(this.min.y, body.position.y - halfHeight);
                this.max.x = Math.max(this.max.x, body.position.x + halfWidth);
                this.max.y = Math.max(this.max.y, body.position.y + halfHeight);
            }
        }
    });

    // Helper Functions ------------------------------------

    /**
     * Flattens the specified array of points onto a unit vector axis,
     * resulting in a one dimensional range of the minimum and
     * maximum value on that axis.
     * @param {Array<game.Vector>} points The points to flatten
     * @param {game.Vector} normal The unit vector axis to flatten on
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
     * @param {game.Vector} aPos The position of the first polygon
     * @param {game.Vector} bPos The position of the second polygon
     * @param {Array<game.Vector>} aPoints The points in the first polygon
     * @param {Array<game.Vector>} bPoints The points in the second polygon
     * @param {game.Vector} axis The axis (unit sized) to test against. The points of both polygons
     *   will be projected onto this axis
     * @param {game.Response=} response A Response object (optional) which will be populated
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
     * @param {game.Vector} line The line segment
     * @param {game.Vector} point The point
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
     * @param {game.Vector} p The point to test
     * @param {game.Circle} c The circle to test
     * @return {Boolean} true if the point is inside the circle, false if it is not
     */
    game.pointInCircle = function pointInCircle(p, c) {
        var differenceV = T_VECTORS.pop().copy(p).subtract(c.position);
        var radiusSq = c.shape.radius * c.shape.radius;
        var distanceSq = differenceV.squaredLength();
        T_VECTORS.push(differenceV);
        // If the distance between is smaller than the radius then the point is inside the circle.
        return distanceSq <= radiusSq;
    }

    /**
     * Check if a point is inside a convex polygon.
     * @param {game.Vector} p The point to test
     * @param {game.Polygon} poly The polygon to test
     * @return {Boolean} true if the point is inside the polygon, false if it is not
     */
    game.pointInPolygon = function pointInPolygon(p, poly) {
        UNIT_SQUARE.position.copy(p);
        T_RESPONSE.clear();
        var result = game.testPolygonPolygon(UNIT_SQUARE, poly, T_RESPONSE);
        if (result) {
            result = T_RESPONSE.aInB;
        }
        return result;
    }

    /**
     * Check if two circles collide.
     * @param {game.Body} a The first circle body
     * @param {game.Body} b The second circle body
     * @param {game.Response=} response Response object (optional) that will be populated if
     *   the circles intersect
     * @return {Boolean} true if the circles intersect, false if they don't
     */
    game.testCircleCircle = function testCircleCircle(a, b, response) {
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
            response.overlapV.copy(differenceV).scale(response.overlap);
            response.aInB = a.shape.radius <= b.shape.radius && dist <= b.shape.radius - a.shape.radius;
            response.bInA = b.shape.radius <= a.shape.radius && dist <= a.shape.radius - b.shape.radius;
        }
        T_VECTORS.push(differenceV);
        return true;
    }

    /**
     * Check if a polygon and a circle collide.
     * @param {game.Polygon} polygon The polygon
     * @param {game.Circle} circle The circle
     * @param {game.Response=} response Response object (optional) that will be populated if
     *   they interset
     * @return {Boolean} true if they intersect, false if they don't
     */
    game.testPolygonCircle = function testPolygonCircle(polygon, circle, response) {
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
            response.overlapV.copy(response.overlapN).scale(response.overlap);
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
     * @param {game.Circle} circle The circle
     * @param {game.Polygon} polygon The polygon
     * @param {game.Response=} response Response object (optional) that will be populated if
     *   they interset
     * @return {Boolean} true if they intersect, false if they don't
     */
    game.testCirclePolygon = function testCirclePolygon(circle, polygon, response) {
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
     * @param {game.Polygon} a The first polygon
     * @param {game.Polygon} b The second polygon
     * @param {game.Response=} response Response object (optional) that will be populated if
     *   they interset
     * @return {Boolean} true if they intersect, false if they don't
     */
    game.testPolygonPolygon = function testPolygonPolygon(a, b, response) {
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
            response.overlapV.copy(response.overlapN).scale(response.overlap);
        }
        return true;
    }

    // Object Pools -----------------------------------------

    /**
     * A pool of `game.Vector` objects that are used in calculations to avoid
     * allocating memory.
     * @type {Array<game.Vector>}
     */
    var T_VECTORS = [];
    for (var i = 0; i < 10; i++) { T_VECTORS.push(new game.Vector()); }

    /**
     * A pool of arrays of numbers used in calculations to avoid allocating
     * memory.
     * @type {Array<Array<Bumber>>}
     */
    var T_ARRAYS = [];
    for (var i = 0; i < 5; i++) { T_ARRAYS.push([]); }

    /**
     * Temporary response used for polygon hit detection.
     * @type {game.Response}
     */
    var T_RESPONSE = new game.Response();

    /**
     * Unit square polygon used for polygon hit detection.
     * @type {game.Polygon}
     */
    var UNIT_SQUARE = new game.Body({
        position: new game.Vector(),
        shape: new game.Rectangle(1, 1).toPolygon()
    });

});
