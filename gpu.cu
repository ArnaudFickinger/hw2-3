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
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

#define NUM_THREADS 256

// Put any static global variables here that you will use throughout the simulation.
int blks;

int* bin_counts_dev;
int* bin_counts_host;

int* prefix_sum_dev;
int* prefix_sum_host;

int* curr_bin_index_dev;
int* curr_bin_index_host;

int* ordered_parts_dev; // each entry is the particle's index in parts
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

__device__ void compute_forces_bin(particle_t* particles, int* prefix_sum_dev, int* ordered_parts_dev, int tid, int bin_num, int num_parts, int num_bins) {

    int curr_offset = prefix_sum_dev[bin_num];
    int stop_offset = (bin_num + 1 >= num_bins) ? num_parts : prefix_sum_dev[bin_num + 1];

    for (int j = curr_offset; j < stop_offset; j++)
        apply_force_gpu(particles[tid], particles[ordered_parts_dev[j]]);
}

__global__ void compute_forces_gpu(particle_t* parts, int num_parts, float size_bin, int num_bins_1d, int* prefix_sum_dev, int* ordered_parts_dev, int num_bins) {
    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts)
        return;

    parts[tid].ax = parts[tid].ay = 0;

    // calculate own bin num
    int bin_x = int(parts[tid].x / size_bin);
    int bin_y = int(parts[tid].y / size_bin);
    int bin_num = bin_x + bin_y * num_bins_1d;

    int left_bin = -1;
    int right_bin = -1;
    int bottom_bin = -1;
    int top_bin = -1;
    int self_bin = bin_num;
    int top_left_bin = -1;
    int top_right_bin = -1;
    int bottom_left_bin = -1;
    int bottom_right_bin = -1;

    if (int(bin_num / num_bins_1d) == 0) { // At top side
        if (bin_num % num_bins_1d == 0) { // top left corner
            right_bin = bin_num + 1;
            bottom_bin = bin_num + num_bins_1d;
            bottom_right_bin = bin_num + num_bins_1d + 1;
            // std::vector<int> vect{ i, i+1, i+num_bins_1d, i + num_bins_1d + 1 };
        } else if (bin_num % num_bins_1d == num_bins_1d - 1) { // top right corner
            left_bin = bin_num - 1;
            bottom_bin = bin_num + num_bins_1d;
            bottom_left_bin = bin_num + num_bins_1d - 1;
            // std::vector<int> vect{ i, i-1, i+num_bins, i + num_bins - 1};
        } else {
            left_bin = bin_num - 1;
            bottom_left_bin = bin_num + num_bins_1d - 1;
            right_bin = bin_num + 1;
            bottom_bin = bin_num + num_bins_1d;
            bottom_right_bin = bin_num + num_bins_1d + 1;
            // std::vector<int> vect{ i - 1, i, i + 1, i - 1 + num_bins, i + num_bins, i + 1 + num_bins };
        }
    } else if (int(bin_num / num_bins_1d) == num_bins_1d - 1) { // bottom side
        if (bin_num % num_bins_1d == 0) { // bottom left corner
            right_bin = bin_num + 1;
            top_bin = bin_num - num_bins_1d;
            top_right_bin = bin_num - num_bins_1d + 1;
            // std::vector<int> vect{ i, i+1, i - num_bins, i - num_bins + 1};
        } else if (bin_num % num_bins_1d == num_bins_1d - 1) { // bottom right corner
            left_bin = bin_num - 1;
            top_bin = bin_num - num_bins_1d;
            top_left_bin = bin_num - num_bins_1d - 1;
            // std::vector<int> vect{ i, i - 1, i - num_bins, i - num_bins - 1};
        } else {
            right_bin = bin_num + 1;
            top_bin = bin_num - num_bins_1d;
            top_right_bin = bin_num - num_bins_1d + 1;
            left_bin = bin_num - 1;
            top_left_bin = bin_num - num_bins_1d - 1;
            // std::vector<int> vect{ i - 1 - num_bins, i - num_bins, i + 1 - num_bins, i - 1, i, i + 1};
        }
    } else {
      if (bin_num % num_bins_1d == 0) { // left side
          right_bin = bin_num + 1;
          top_bin = bin_num - num_bins_1d;
          top_right_bin = bin_num - num_bins_1d + 1;
          bottom_bin = bin_num + num_bins_1d;
          bottom_right_bin = bin_num + num_bins_1d + 1;
          // std::vector<int> vect{ i - num_bins, i+1 - num_bins, i, i + 1, i + num_bins, i + 1 + num_bins};
      } else if (bin_num % num_bins_1d == num_bins_1d - 1) { // right side
          left_bin = bin_num - 1;
          top_bin = bin_num - num_bins_1d;
          top_left_bin = bin_num - num_bins_1d - 1;
          bottom_bin = bin_num + num_bins_1d;
          bottom_left_bin = bin_num + num_bins_1d - 1;
          // std::vector<int> vect{ i - 1 - num_bins, i - num_bins, i, i - 1, i - 1 + num_bins, i + num_bins};
      } else {
          right_bin = bin_num + 1;
          top_bin = bin_num - num_bins_1d;
          top_right_bin = bin_num - num_bins_1d + 1;
          bottom_bin = bin_num + num_bins_1d;
          bottom_right_bin = bin_num + num_bins_1d + 1;
          left_bin = bin_num - 1;
          top_left_bin = bin_num - num_bins_1d - 1;
          bottom_left_bin = bin_num + num_bins_1d - 1;
          // std::vector<int> vect{ i - 1 - num_bins, i - num_bins, i + 1 - num_bins, i - 1, i, i + 1, i - 1 + num_bins, i + num_bins, i + 1 + num_bins};
      }
    }

    if (left_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, left_bin, num_parts, num_bins);
    }

    if (right_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, right_bin, num_parts, num_bins);
    }

    if (bottom_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, bottom_bin, num_parts, num_bins);
    }

    if (top_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, top_bin, num_parts, num_bins);
    }

    compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, self_bin, num_parts, num_bins);

    if (top_left_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, top_left_bin, num_parts, num_bins);
    }

    if (top_right_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, top_right_bin, num_parts, num_bins);
    }

    if (bottom_left_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, bottom_left_bin, num_parts, num_bins);
    }

    if (bottom_right_bin != -1) {
        compute_forces_bin(parts, prefix_sum_dev, ordered_parts_dev, tid, bottom_right_bin, num_parts, num_bins);
    }



    // find up to 9 relevant bins
    // compute_forces_bin for each neighboring bin

    // for (int j = 0; j < num_parts; j++)
    //     apply_force_gpu(particles[tid], particles[j]);
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

    // int* bin_counts_dev;
    // int* bin_counts_host;
    bin_counts_host = (int*) calloc(num_bins, sizeof(int));
    cudaMalloc((void**) &bin_counts_dev, sizeof(int) * num_bins);
    cudaMemset(bin_counts_dev, 0, num_bins * sizeof(int));

    // int* prefix_sum_dev;
    // int* prefix_sum_host;
    prefix_sum_host = (int*) calloc(num_bins+1, sizeof(int));
    cudaMalloc((void**) &prefix_sum_dev, sizeof(int) * (num_bins + 1));

    // int* curr_bin_index_dev;
    // int* curr_bin_index_host;
    curr_bin_index_host = (int*) calloc(num_bins+1, sizeof(int));
    cudaMalloc((void**) &curr_bin_index_dev, sizeof(int) * (num_bins + 1));

    // int* ordered_parts_dev;
    // int* ordered_parts_host;
    ordered_parts_host = (int*) calloc(num_parts, sizeof(int));
    cudaMalloc((void**) &ordered_parts_dev, sizeof(int) * num_parts);
}

__global__ void update_bin_counts(particle_t* parts, int num_parts, int* bin_counts, float size_bin, int num_bins_1d, int num_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts) {
      return;
    }

    int bin_x = int(parts[tid].x / size_bin);
    int bin_y = int(parts[tid].y / size_bin);
    int bin_num = bin_x + bin_y * num_bins_1d;

    atomicAdd(&bin_counts[bin_num], 1);
}

__global__ void order_particles(particle_t* parts, int num_parts, float size_bin, int num_bins_1d, int* curr_bin_index_dev, int* ordered_parts_dev) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts) {
      return;
    }

    int bin_x = int(parts[tid].x / size_bin);
    int bin_y = int(parts[tid].y / size_bin);
    int bin_num = bin_x + bin_y * num_bins_1d;

    int index = atomicAdd(&curr_bin_index_dev[bin_num], 1);
    ordered_parts_dev[index] = tid;
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // parts live in GPU memory
    // Rewrite this function

    update_bin_counts<<<blks, NUM_THREADS>>>(parts, num_parts, bin_counts_dev, size_bin, num_bins_1d, num_bins);

    // DEBUG BIN COUNTS
    // std::cout << "BIN COUNTS" << std::endl;
    // cudaMemcpy(bin_counts_host, bin_counts_dev, sizeof(int) * num_bins, cudaMemcpyDeviceToHost);
    // for (int i = 0; i < num_bins; i++) {
    //     std::cout << bin_counts_host[i] << std::endl;
    // }

    // cudaMemcpy(prefix_sum_dev, bin_counts_dev, num_bins * sizeof(int), cudaMemcpyDeviceToDevice);

    thrust::exclusive_scan(thrust::device, bin_counts_dev, bin_counts_dev + num_bins, prefix_sum_dev, 0);

    // DEBUG PREFIX SUM
    // std::cout << "PREFIX SUM" << std::endl;
    // cudaMemcpy(prefix_sum_host, prefix_sum_dev, sizeof(int) * (num_bins + 1), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < num_bins + 1; i++) {
    //     std::cout << prefix_sum_host[i] << std::endl;
    // }

    cudaMemcpy(curr_bin_index_dev, prefix_sum_dev, (num_bins + 1) * sizeof(int), cudaMemcpyDeviceToDevice);

    // DEBUG PRE-ORDERED CURR INDEXES
    // std::cout << "PRE-ORDERED CURR INDEXES" << std::endl;
    // cudaMemcpy(curr_bin_index_host, curr_bin_index_dev, sizeof(int) * (num_bins + 1), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < num_bins + 1; i++) {
    //     std::cout << curr_bin_index_host[i] << std::endl;
    // }

    order_particles<<<blks, NUM_THREADS>>>(parts, num_parts, size_bin, num_bins_1d, curr_bin_index_dev, ordered_parts_dev);

    // DEBUG ORDERED PARTS INDEXES
    // std::cout << "ORDERED PARTS" << std::endl;
    // cudaMemcpy(ordered_parts_host, ordered_parts_dev, sizeof(int) * num_parts, cudaMemcpyDeviceToHost);
    // for (int i = 0; i < num_parts; i++) {
    //     std::cout << ordered_parts_host[i] << std::endl;
    // }
    //
    // std::cout << "ORDERED CURR INDEXES" << std::endl;
    // cudaMemcpy(curr_bin_index_host, curr_bin_index_dev, sizeof(int) * (num_bins + 1), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < num_bins + 1; i++) {
    //     std::cout << curr_bin_index_host[i] << std::endl;
    // }

    // Compute forces
    compute_forces_gpu<<<blks, NUM_THREADS>>>(parts, num_parts, size_bin, num_bins_1d, prefix_sum_dev, ordered_parts_dev, num_bins);

    // Move particles
    move_gpu<<<blks, NUM_THREADS>>>(parts, num_parts, size);

    cudaMemset(bin_counts_dev, 0, num_bins * sizeof(int));
    // std::cout << "end step" << std::endl;
}
