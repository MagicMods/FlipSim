import { eventBus } from '../../util/eventManager.js';

export class CollisionSystem {
  constructor({
    enabled = true,
    gridSize = 10,
    repulsion = 0.5,
    damping = 0.98,
    particleRestitution = 0.8,
    particleRadius = 0.01,
  } = {}) {
    this.enabled = enabled;
    this.gridSize = gridSize;
    this.cellSize = 1 / gridSize;
    this.repulsion = repulsion;
    this.damping = damping;
    this.particleRestitution = particleRestitution;
    this.particleRadius = particleRadius;
    this.maxParticleRadius = particleRadius;

    // INITIALIZE THIS - it was missing
    this.grid = [];
    for (let i = 0; i < gridSize * gridSize; i++) {
      this.grid.push([]);
    }

    this.particleSystem = null;

    // Assuming initializeGrid exists and needs to be called
    this.initializeGrid(); // Call initial grid setup

    // Subscribe to parameter updates
    eventBus.on('simParamsUpdated', this.handleParamsUpdate.bind(this));
  }

  // Add handler for simParams updates
  handleParamsUpdate({ simParams }) {
    // Update base radius first
    if (simParams?.simulation?.particleRadius !== undefined) {
      this.particleRadius = simParams.simulation.particleRadius;
    }

    // Update collision parameters
    if (simParams?.collision) {
      const previousGridSize = this.gridSize; // Store previous size for comparison

      this.enabled = simParams.collision.enabled ?? this.enabled;
      this.gridSize = simParams.collision.gridSize ?? this.gridSize;
      this.repulsion = simParams.collision.repulsion ?? this.repulsion;
      this.particleRestitution = simParams.collision.particleRestitution ?? this.particleRestitution;

      // Re-initialize grid if gridSize has changed
      if (this.gridSize !== previousGridSize && typeof this.initializeGrid === 'function') {
        console.log(`CollisionSystem: Grid size changed from ${previousGridSize} to ${this.gridSize}. Reinitializing grid.`);
        this.initializeGrid(); // Assuming this method exists to handle grid resizing
      }
    }

    // Update maxParticleRadius based on turbulence settings AFTER updating base radius
    if (simParams?.turbulence?.affectScale === true) {
      // If affectScale is true, use the maxScale from turbulence params, fallback to current base radius
      this.maxParticleRadius = simParams.turbulence.maxScale ?? this.particleRadius;
    } else {
      // If affectScale is false or turbulence params are missing, max radius is just the base radius
      this.maxParticleRadius = this.particleRadius;
    }
    // console.log(`CollisionSystem updated: enabled=${this.enabled}, gridSize=${this.gridSize}, baseRadius=${this.particleRadius}, maxRadius=${this.maxParticleRadius}`);
  }

  // Add initializeGrid method if it doesn't exist (or update existing one)
  initializeGrid() {
    this.cellSize = 1 / this.gridSize;
    // Reset grid array based on new size
    this.grid = [];
    const gridCellCount = this.gridSize * this.gridSize;
    for (let i = 0; i < gridCellCount; i++) {
      this.grid.push([]);
    }
    // console.log(`Collision grid initialized with size ${this.gridSize}x${this.gridSize}`); // REMOVED
  }

  update(particles, velocitiesX, velocitiesY) {
    if (!this.enabled) return;

    this.updateGrid(particles);
    this.resolveCollisions(particles, velocitiesX, velocitiesY);
  }

  updateGrid(particles) {
    // Reset grid at the start of each update
    this.grid = [];
    for (let i = 0; i < this.gridSize * this.gridSize; i++) {
      this.grid.push([]);
    }

    // Then add particles to grid cells
    const numParticles = particles.length / 2;
    for (let i = 0; i < numParticles; i++) {
      const x = particles[i * 2];
      const y = particles[i * 2 + 1];

      // Ensure coordinates are within bounds
      if (x < 0 || x > 1 || y < 0 || y > 1) continue;

      const cellX = Math.floor(x * this.gridSize);
      const cellY = Math.floor(y * this.gridSize);
      const cellIndex = cellY * this.gridSize + cellX;

      // Ensure the cell exists before pushing
      if (this.grid[cellIndex]) {
        this.grid[cellIndex].push(i);
      }
    }
  }

  resolveCollisions(particles, velocitiesX, velocitiesY) {
    for (let cellIndex = 0; cellIndex < this.grid.length; cellIndex++) {
      this.checkCellCollisions(cellIndex, particles, velocitiesX, velocitiesY);
    }
  }

  // --- Add helper method to calculate search radius --- 
  calculateSearchRadius() {
    if (!this.cellSize || this.cellSize <= 0) return 1; // Safety check
    // Ensure maxParticleRadius has a fallback
    const maxRadius = this.maxParticleRadius || this.particleRadius || 0.01;
    // Calculate search radius in cells, minimum 1
    return Math.max(1, Math.ceil(maxRadius / this.cellSize));
  }
  // --- End helper method ---

  checkCellCollisions(cellIndex, particles, velocitiesX, velocitiesY) {
    const cell = this.grid[cellIndex];
    if (!cell) return; // Add safety check for empty/invalid cell

    const x = cellIndex % this.gridSize;
    const y = Math.floor(cellIndex / this.gridSize);
    const searchRadius = this.calculateSearchRadius(); // Calculate dynamic search radius

    for (let i = 0; i < cell.length; i++) {
      const particleI = cell[i];

      // Same cell collision check (inner loop)
      for (let j = i + 1; j < cell.length; j++) {
        this.resolveCollision(
          particleI,
          cell[j],
          particles,
          velocitiesX,
          velocitiesY
        );
      }

      // Check neighboring cells using nested loops and searchRadius
      for (let offsetY = -searchRadius; offsetY <= searchRadius; offsetY++) {
        for (let offsetX = -searchRadius; offsetX <= searchRadius; offsetX++) {
          // Skip checking the cell against itself OR checking previously checked pairs (e.g., only check positive offsets relative to current)
          // Simplification: Only check offsets where interaction hasn't been checked by the 'other' cell
          // This avoids double checks. We check x+1,y; x,y+1; x+1,y+1; x-1,y+1; x+1, y-1 (relative positive offsets)
          // A more robust way is needed, but let's keep it simple first
          // For now, we'll check all neighbors in the square and rely on resolveCollision check
          if (offsetX === 0 && offsetY === 0) continue; // Skip self

          const neighborX = x + offsetX;
          const neighborY = y + offsetY;

          // Check bounds before calling checkNeighborCell
          if (neighborX >= 0 && neighborX < this.gridSize && neighborY >= 0 && neighborY < this.gridSize) {
            this.checkNeighborCell(
              neighborX,
              neighborY,
              particleI,
              particles,
              velocitiesX,
              velocitiesY
            );
          }
        }
      }
    }
  }

  checkNeighborCell(x, y, particleI, particles, velocitiesX, velocitiesY) {
    if (x < 0 || x >= this.gridSize || y < 0 || y >= this.gridSize) {
      return;
    }

    const neighborIndex = y * this.gridSize + x;
    const neighborCell = this.grid[neighborIndex];

    for (const particleJ of neighborCell) {
      this.resolveCollision(
        particleI,
        particleJ,
        particles,
        velocitiesX,
        velocitiesY
      );
    }
  }

  resolveCollision(i, j, particles, velocitiesX, velocitiesY) {
    const dx = particles[j * 2] - particles[i * 2];
    const dy = particles[j * 2 + 1] - particles[i * 2 + 1];
    const distSq = dx * dx + dy * dy;

    // Replace existing radius logic with direct access to particleRadii array
    let radiusI = (this.particleSystem?.particleRadii?.[i]) ?? this.particleRadius;
    let radiusJ = (this.particleSystem?.particleRadii?.[j]) ?? this.particleRadius;

    // Access rest density from particleSystem reference
    let densityFactor = 1.0;
    if (this.particleSystem && this.particleSystem.restDensity) {
      densityFactor = Math.max(0.5, Math.min(2.0, this.particleSystem.restDensity / 3));
    }

    // Modify minimum separation distance based on rest density
    // Higher density = particles should be closer (smaller minDist)
    const minDist = (radiusI + radiusJ) * (1 / densityFactor);

    const minDistSq = minDist * minDist;

    if (distSq < minDistSq) {
      const dist = Math.sqrt(distSq);
      const nx = dx / dist;
      const ny = dy / dist;

      // Calculate overlap
      const overlap = minDist - dist;

      // Convert repulsion into a smooth force
      const repulsionForce = Math.pow(overlap / minDist, 1.5) * this.repulsion;

      // Apply repulsion as velocity change instead of position correction
      const repulsionX = nx * repulsionForce;
      const repulsionY = ny * repulsionForce;

      // Apply smoother velocity changes
      velocitiesX[i] -= repulsionX;
      velocitiesY[i] -= repulsionY;
      velocitiesX[j] += repulsionX;
      velocitiesY[j] += repulsionY;

      // Handle collision response only for significant relative velocity
      const vx = velocitiesX[j] - velocitiesX[i];
      const vy = velocitiesY[j] - velocitiesY[i];
      const vn = vx * nx + vy * ny;

      if (vn < 0) {
        // Softer restitution based on overlap
        const restitution =
          this.particleRestitution * (1.0 - overlap / minDist);
        const impulse = -(1 + restitution) * vn * 0.5;

        velocitiesX[i] -= nx * impulse;
        velocitiesY[i] -= ny * impulse;
        velocitiesX[j] += nx * impulse;
        velocitiesY[j] += ny * impulse;
      }

      // Apply minimal position correction to prevent extreme overlap
      const correction = Math.max(0, overlap - minDist * 0.1) * 0.5;
      particles[i * 2] -= nx * correction;
      particles[i * 2 + 1] -= ny * correction;
      particles[j * 2] += nx * correction;
      particles[j * 2 + 1] += ny * correction;
    }
  }

  applyImpulse(i, j, nx, ny, vn, velocitiesX, velocitiesY) {
    const impulse = -(1 + this.particleRestitution) * vn * 0.5;

    velocitiesX[i] -= impulse * nx;
    velocitiesY[i] -= impulse * ny;
    velocitiesX[j] += impulse * nx;
    velocitiesY[j] += impulse * ny;

    // Apply damping
    velocitiesX[i] *= this.damping;
    velocitiesY[i] *= this.damping;
    velocitiesX[j] *= this.damping;
    velocitiesY[j] *= this.damping;
  }

  applyRepulsion(i, j, nx, ny, dist, velocitiesX, velocitiesY) {
    const overlap = this.particleRadius * 2 - dist;
    const repulsionForce = overlap * this.repulsion;

    velocitiesX[i] -= nx * repulsionForce;
    velocitiesY[i] -= ny * repulsionForce;
    velocitiesX[j] += nx * repulsionForce;
    velocitiesY[j] += ny * repulsionForce;
  }

  reset() {
    // Clear all grid cells
    this.grid.forEach((cell) => (cell.length = 0));
  }
}