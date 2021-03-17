#include "common.h"
#include <cuda.h>

#define NUM_THREADS 256

// Put any static global variables here that you will use throughout the simulation.
int blks;

__device__ void apply_force_gpu(particle_t& particle, particle_t& neighbor) {
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
        return;
    // r2 = fmax( r2, min_r*min_r );
    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r;
    double r = sqrt(r2);

    //
    //  very simple short-range repulsive force
    //
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

__global__ void compute_forces_gpu(particle_t* particles, int num_parts) {
    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts)
        return;

    particles[tid].ax = particles[tid].ay = 0;
    for (int j = 0; j < num_parts; j++)
        apply_force_gpu(particles[tid], particles[j]);
}

__global__ void move_gpu(particle_t* particles, int num_parts, double size) {

    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts)
        return;

    particle_t* p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x += p->vx * dt;
    p->y += p->vy * dt;

    //
    //  bounce from walls
    //
    while (p->x < 0 || p->x > size) {
        p->x = p->x < 0 ? -(p->x) : 2 * size - p->x;
        p->vx = -(p->vx);
    }
    while (p->y < 0 || p->y > size) {
        p->y = p->y < 0 ? -(p->y) : 2 * size - p->y;
        p->vy = -(p->vy);
    }
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    // You can use this space to initialize data objects that you may need
    // This function will be called once before the algorithm begins
    // parts live in GPU memory
    // Do not do any particle simulation here

    blks = (num_parts + NUM_THREADS - 1) / NUM_THREADS;

    /////
    num_bins_side = int(size / cutoff);
    num_bins = num_bins_side * num_bins_side;
    size_bin = size / num_bins_side;
    /////
}

__global__ void create_bin_counts(particle_t* parts, int num_parts, int* bin_counts, int size_bin, int num_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts) {
      return;
    }
    // particle_t& part = parts[tid];

    // int bin_x = int(part.x / size_bin);
    // int bin_y = int(part.y / size_bin);
    // int bin_num = bin_x + bin_y * num_bins;

    int bin_x = int(parts[tid].x / size_bin);
    int bin_y = int(parts[tid].y / size_bin);
    int bin_num = bin_x + bin_y * num_bins;

    // int* cpu_bin_num = (int*) malloc(sizeof(int));
    // cudaMemcpy(cpu_bin_num, bin_num, sizeof(int), cudaMemcpyDeviceToHost);

    // std::cout << cpu_bin_num << ",\t";

    bin_counts[bin_num] = 5;

    // atomicAdd(&bin_counts[bin_num], 1);
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // parts live in GPU memory
    // Rewrite this function

    /////
    if (!parts_ordered_inds) {
        cudaMalloc((void**) &parts_ordered_inds, size * sizeof(int));
    }

    if (!bin_counts) {
        cudaMalloc((void**) &bin_counts, num_bins * sizeof(int));
    }

    cudaMemset(parts_ordered_inds, -1, size * sizeof(int));
    cudaMemset(bin_counts, 0, num_bins * sizeof(int));

    int* cpu_bin_counts = (int*) malloc(num_bins * sizeof(int));
    cudaMemcpy(cpu_bin_counts, bin_counts, num_bins * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < 2; i++) {
        std::cout << i << ": "<<  cpu_bin_counts[i] << ",\t";
    }

    create_bin_counts<<<blks, NUM_THREADS>>>(parts, num_parts, bin_counts, size_bin, num_bins);
    /////

    // Compute forces
    compute_forces_gpu<<<blks, NUM_THREADS>>>(parts, num_parts);

    // Move particles
    move_gpu<<<blks, NUM_THREADS>>>(parts, num_parts, size);
}
