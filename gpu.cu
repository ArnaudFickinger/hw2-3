// #include "common.h"
// #include <cuda.h>
// #include <iostream>
// #include <thrust/scan.h>
// #include <thrust/execution_policy.h>
//
// #include <thrust/scan.h>
// #include <thrust/execution_policy.h>
//
// #define NUM_THREADS 256
//
// // Put any static global variables here that you will use throughout the simulation.
// int blks;
// int num_bins_1d;
// int num_bins;
// float size_bin;
// int size_bin_counts;
//
// int* bin_counts_host;
// int* bin_counts_host_check;
// int* bin_counts_device;
//
// int* bin_counts_sum_device;
//
// int* bin_counts_incremental_host;
// __device__ int* bin_counts_incremental_device;
//
// int* ordered_particles_host;
// __device__ int* ordered_particles_device;
//
// int* bin_counts_sum_check;
//
// int sim_number = 0;
//
// __device__ void apply_force_gpu(particle_t& particle, particle_t& neighbor) {
//     double dx = neighbor.x - particle.x;
//     double dy = neighbor.y - particle.y;
//     double r2 = dx * dx + dy * dy;
//     if (r2 > cutoff * cutoff)
//         return;
//     // r2 = fmax( r2, min_r*min_r );
//     r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r;
//     double r = sqrt(r2);
//
//     //
//     //  very simple short-range repulsive force
//     //
//     double coef = (1 - cutoff / r) / r2 / mass;
//     particle.ax += coef * dx;
//     particle.ay += coef * dy;
// }
//
// __global__ void compute_forces_gpu_naive(particle_t* particles, int num_parts) {
//     // Get thread (particle) ID
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid >= num_parts)
//         return;
//
//     particles[tid].ax = particles[tid].ay = 0;
//     for (int j = 0; j < num_parts; j++)
//         apply_force_gpu(particles[tid], particles[j]);
// }
//
// __global__ void compute_forces_gpu(particle_t* ordered_particles, int* bin_counts_sum, int num_parts_, float size_bin_, int num_bins_1d_, int num_bins_) {
//     // Get thread (particle) ID
//     bin_counts_sum[num_bins_] = num_parts_;
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid >= num_parts_)
//         return;
//
//     int bin_x = int(ordered_particles[tid].x / size_bin_);
//     int bin_y = int(ordered_particles[tid].y / size_bin_);
//     int bin_num = bin_x + bin_y * num_bins_1d_;
//
//     int index_first = 0;
//     int index_last = 0;
//
//     if (bin_num==0)
//        index_first = 0;
//     else
//        index_first = bin_counts_sum[bin_num-1];
//     if (bin_num>=num_bins_-2)
//        index_last = num_parts_;
//     else
//        index_last = bin_counts_sum[bin_num+2];
//
//     ordered_particles[tid].ax = ordered_particles[tid].ay = 0;
//     for (int j = index_first; j < index_last; j++)
//         apply_force_gpu(ordered_particles[tid], ordered_particles[j]);
// }
//
// __global__ void move_gpu(particle_t* particles, int num_parts, double size) {
//
//     // Get thread (particle) ID
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid >= num_parts)
//         return;
//
//     particle_t* p = &particles[tid];
//     //
//     //  slightly simplified Velocity Verlet integration
//     //  conserves energy better than explicit Euler method
//     //
//     p->vx += p->ax * dt;
//     p->vy += p->ay * dt;
//     p->x += p->vx * dt;
//     p->y += p->vy * dt;
//
//     //
//     //  bounce from walls
//     //
//     while (p->x < 0 || p->x > size) {
//         p->x = p->x < 0 ? -(p->x) : 2 * size - p->x;
//         p->vx = -(p->vx);
//     }
//     while (p->y < 0 || p->y > size) {
//         p->y = p->y < 0 ? -(p->y) : 2 * size - p->y;
//         p->vy = -(p->vy);
//     }
// }
//
// __global__ void update_bin_counts(particle_t* parts, int num_parts, int* bin_counts, float size_bin_, int num_bins_) {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid >= num_parts) {
//       return;
//     }
//
//     int bin_x = int(parts[tid].x / size_bin_);
//     int bin_y = int(parts[tid].y / size_bin_);
//     int bin_num = bin_x + bin_y * num_bins_;
//
//     atomicAdd(&bin_counts[bin_num], 1);
// }
//
// __global__ void order_particle(particle_t* parts, int num_parts, float size_bin_, int num_bins_, int* bin_counts_incremental_device_, int* ordered_particles_device_) {
//
//
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid >= num_parts) {
//       return;
//     }
//
//     int bin_x = int(parts[tid].x / size_bin_);
//     int bin_y = int(parts[tid].y / size_bin_);
//     int bin_num = bin_x + bin_y * num_bins_;
//
//     // needs to be atomic
//     // __syncthreads();
//     // int index = bin_counts_incremental_device_[bin_num] - 1;
//     // cudaMemset(&ordered_particles_device_[index], tid, sizeof(int));
//
//     int index = bin_counts_incremental_device_[bin_num];
//     ordered_particles_device_[index] = tid;
//     atomicAdd(&bin_counts_incremental_device_[bin_num], 1);
//
//     // int index = bin_counts_incremental_device_[bin_num];
//     // int i = 0;
//     // while (true) {
//     //     // if (ordered_particles_device_[index] == -1) {
//     //         ordered_particles_device_[index] = tid;
//     //         // break;
//     //     // } else {
//     //         // index++;
//     //     // }
//     //     if (i == 10) {
//     //         break;
//     //     }
//     //     i++;
//     // }
//
// }
//
// // __global__ void update_bin_counts_test(particle_t *parts, int num_parts, int *bin_counts, float size_bin_, int num_bins_) {
// //
// //     // bin_counts[0]=1;
// //     atomicAdd(&bin_counts[0], 1);
// //     // bin_counts[bin_num]+=1;
// //
// //     // atomicAdd(&bin_counts[bin_num], 1);
// // }
//
// void init_simulation(particle_t* parts, int num_parts, double size) {
//     // You can use this space to initialize data objects that you may need
//     // This function will be called once before the algorithm begins
//     // parts live in GPU memory
//     // Do not do any particle simulation here
//
//     blks = (num_parts + NUM_THREADS - 1) / NUM_THREADS;
//
//     num_bins_1d = 2;
//     size_bin = size/num_bins_1d;
//     num_bins = num_bins_1d*num_bins_1d;
//     size_bin_counts = num_bins* sizeof(int);
//     bin_counts_host = (int*)calloc(num_bins, sizeof(int));
//     bin_counts_host_check = (int*)calloc(num_bins, sizeof(int));
//
//
//     cudaMalloc(&bin_counts_device, size_bin_counts);
//     cudaMemcpy(bin_counts_device, bin_counts_host, size_bin_counts, cudaMemcpyHostToDevice);
//
//
//     cudaMalloc((void**)&ordered_particles_host, num_parts * sizeof(int));
//     cudaMemcpyToSymbol(ordered_particles_device, &ordered_particles_host, sizeof(int*));
//     cudaMemset(ordered_particles_host, -1, num_parts * sizeof(int));
//
//     // int* bin_counts_incremental_host;
//     // __device__ int* bin_counts_incremental_device;
//
//     cudaMalloc((void**) &bin_counts_sum_device, (num_bins + 1) * sizeof(int));
//     // cudaMalloc((void**) &bin_counts_incremental_device, (num_bins + 1) * sizeof(int));
//     cudaMalloc((void**) &bin_counts_incremental_host, (num_bins + 1) * sizeof(int));
//     cudaMemcpyToSymbol(bin_counts_incremental_device, &bin_counts_incremental_host, sizeof(int*));
//
//
//     bin_counts_sum_check = (int*) malloc(sizeof(int) * (num_bins + 1));
//
// }
//
//
// void simulate_one_step(particle_t* parts, int num_parts, double size) {
//     // parts live in GPU memory
//     // Rewrite this function
//
//     cudaMemset(bin_counts_device, 0, num_bins * sizeof(int));
//
//     update_bin_counts<<<blks, NUM_THREADS>>>(parts, num_parts, bin_counts_device, size_bin, num_bins_1d);
//
//     // cudaMemcpy(bin_counts_host_check, bin_counts_device, size_bin_counts, cudaMemcpyDeviceToHost);
//     //
//     // std::cout << "new step" << ",\t";
//     //
//     // for (int i = 0; i < num_bins; i++) {
//     //     std::cout << i << ": "<<  bin_counts_host_check[i] << ",\t";
//     // }
//
//     thrust::exclusive_scan(thrust::device, bin_counts_device, bin_counts_device + num_bins, bin_counts_sum_device, 0);
//
//     cudaMemcpy(bin_counts_sum_check, bin_counts_sum_device, (num_bins + 1) * sizeof(int), cudaMemcpyDeviceToHost);
//
//     std::cout << "new step" << ",\t";
//
//     for (int i = 0; i < num_bins + 1; i++) {
//         std::cout << i << ": "<<  bin_counts_sum_check[i] << ",\t";
//     }
//
//     order_particle<<<blks, NUM_THREADS>>>(parts, num_parts, size_bin, num_bins, bin_counts_incremental_host, ordered_particles_host);
//
//     if (sim_number == 0) {
//
//         int* ordered_parts_check = (int*) malloc(sizeof(int) * num_parts);
//
//         cudaMemcpy(ordered_parts_check, ordered_particles_device, (num_parts) * sizeof(int), cudaMemcpyDeviceToHost);
//
//         for (int i = 0; i < num_parts; i++) {
//             std::cout << i << ": "<<  ordered_parts_check[i] << ",\t";
//         }
//     }
//
//
//     sim_number++;
//     // sort array
//     // cudaMemcpy(bin_counts_incremental_device, bin_counts_sum_device, size_bin_counts, cudaMemcpyDeviceToDevice);
//     // order_particle<<blks, NUM_THREADS>>(parts, num_parts, bin_counts_device, size_bin, num_bins, bin_counts_incremental_device, ordered_particles_device);
//     //
//
//     // // Compute forces
//     // compute_forces_gpu<<<blks, NUM_THREADS>>>(parts, num_parts);
//     //
//     // // Move particles
//     // move_gpu<<<blks, NUM_THREADS>>>(parts, num_parts, size);
// }









#include "common.h"
#include <cuda.h>
#include <iostream>

#define NUM_THREADS 256

// Put any static global variables here that you will use throughout the simulation.
int blks;

__device__ int* bin_counts_dev;
int* bin_counts_host;

__device__ int* prefix_sum_dev;
int* prefix_sum_host;

__device__ int* curr_bin_index_dev;
int* curr_bin_index_host;

__device__ int* ordered_parts_dev; // each entry is the particle's index in parts
int* ordered_parts_host;

int num_bins_1d;
int num_bins;
float size_bin;

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

    num_bins_1d = 2;
    size_bin = size / num_bins_1d;
    num_bins = num_bins_1d * num_bins_1d;

    // __device__ int* bin_counts_dev;
    // int* bin_counts_host;
    cudaMalloc((void**) &bin_counts_host, num_bins * sizeof(int));
    cudaMemcpyToSymbol(bin_counts_dev, &bin_counts_host, sizeof(int *));
    cudaMemset(bin_counts_host, 0, num_bins * sizeof(int));

    // __device__ int* prefix_sum_dev;
    // int* prefix_sum_host;
    cudaMalloc((void**) &prefix_sum_host, (num_bins + 1) * sizeof(int));
    cudaMemcpyToSymbol(prefix_sum_dev, &prefix_sum_host, sizeof(int *));
    cudaMemset(prefix_sum_host, 0, (num_bins + 1) * sizeof(int));

    // __device__ int* curr_bin_index_dev;
    // int* curr_bin_index_host;
    cudaMalloc((void**) &curr_bin_index_host, (num_bins + 1) * sizeof(int));
    cudaMemcpyToSymbol(curr_bin_index_dev, &curr_bin_index_host, sizeof(int *));
    cudaMemset(curr_bin_index_host, 0, (num_bins + 1) * sizeof(int));

    // __device__ int* ordered_parts_dev;
    // int* ordered_parts_host;
    cudaMalloc((void**) &ordered_parts_host, num_parts * sizeof(int));
    cudaMemcpyToSymbol(ordered_parts_dev, &ordered_parts_host, sizeof(int *));
    cudaMemset(ordered_parts_host, 0, num_parts * sizeof(int));
}

__global__ void update_bin_counts(particle_t* parts, int num_parts, int* bin_counts, float size_bin, int num_bins_1d, int num_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts) {
      return;
    }

    int bin_x = int(parts[tid].x / size_bin);
    int bin_y = int(parts[tid].y / size_bin);
    int bin_num = bin_x + bin_y * num_bins;

    atomicAdd(&bin_counts[bin_num], 1);
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // parts live in GPU memory
    // Rewrite this function

    update_bin_counts<<<blks, NUM_THREADS>>>(parts, num_parts, bin_counts_host, size_bin, num_bins_1d, num_bins);

    std::cout << bin_counts_host[0] << std::endl;
    // for (int i = 0; i < num_bins; i++) {
    //     std::cout << bin_counts_host[i] << std::endl;
    // }

    // Compute forces
    // compute_forces_gpu<<<blks, NUM_THREADS>>>(parts, num_parts);

    // Move particles
    // move_gpu<<<blks, NUM_THREADS>>>(parts, num_parts, size);
}
