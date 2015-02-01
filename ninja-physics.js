game.module(
    'plugins.ninja-physics'
)
.require(
    'engine.physics'
)
.body(function() { 'use strict';

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

});
