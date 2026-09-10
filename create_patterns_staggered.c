#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(path) _mkdir(path)
#else
#define MAKE_DIR(path) mkdir(path, 0777)
#endif

typedef struct {
    int r1, r2, c1, c2; // Row min/max, Col min/max in sub-pixel units
} Obstacle;

// Render obstacles to discrete integer sub-pixel grid
void render_to_grid(Obstacle *obstacles, int obs_count, int *grid, int grid_dim) {
    memset(grid, 0, grid_dim * grid_dim * sizeof(int));

    for (int k = 0; k < obs_count; k++) {
        for (int r = obstacles[k].r1; r <= obstacles[k].r2; r++) {
            for (int c = obstacles[k].c1; c <= obstacles[k].c2; c++) {
                if (r >= 0 && r < grid_dim && c >= 0 && c < grid_dim) {
                    grid[r * grid_dim + c] = 1;
                }
            }
        }
    }
}

// Export TXT and PGM files for the scaled grid
void save_pattern(int *grid, int grid_dim, const char *dir_name, const char *pattern_name) {
    char txt_path[256], pgm_path[256];
    sprintf(txt_path, "%s/%s.txt", dir_name, pattern_name);
    sprintf(pgm_path, "%s/%s.pgm", dir_name, pattern_name);

    // Save TXT (1-based occupied pixels)
    FILE *txt_f = fopen(txt_path, "w");
    if (txt_f) {
        int pixel_count = 0;
        for (int i = 0; i < grid_dim * grid_dim; i++) {
            if (grid[i]) pixel_count++;
        }
        fprintf(txt_f, "%d\n", pixel_count);
        for (int r = 0; r < grid_dim; r++) {
            for (int c = 0; c < grid_dim; c++) {
                if (grid[r * grid_dim + c]) {
                    fprintf(txt_f, "%d %d %d %d\n", r + 1, r + 1, c + 1, c + 1);
                }
            }
        }
        fclose(txt_f);
    }

    // Save PGM Image
    FILE *pgm_f = fopen(pgm_path, "w");
    if (pgm_f) {
        fprintf(pgm_f, "P2\n%d %d\n255\n", grid_dim, grid_dim);
        for (int r = 0; r < grid_dim; r++) {
            for (int c = 0; c < grid_dim; c++) {
                fprintf(pgm_f, "%d ", grid[r * grid_dim + c] ? 0 : 255);
            }
            fprintf(pgm_f, "\n");
        }
        fclose(pgm_f);
    }
}

int main() {
    int N_base, shift_choice, num_random_patterns;
    int base_spacing_x, base_spacing_y;

    printf("Enter base N_grid (e.g., 40): ");
    if (scanf("%d", &N_base) != 1) return 1;

    printf("Select shift type:\n 1) Quarter pixel (0.25)\n 2) Half pixel (0.50)\n 3) Full pixel (1.00)\nChoice: ");
    if (scanf("%d", &shift_choice) != 1) return 1;

    int scale_factor;
    if (shift_choice == 1)      scale_factor = 4; // 1 unit shift = 1/4 pixel
    else if (shift_choice == 2) scale_factor = 2; // 1 unit shift = 1/2 pixel
    else                        scale_factor = 1; // 1 unit shift = 1 full pixel

    printf("Enter spacing X in base units (e.g., 2): ");
    if (scanf("%d", &base_spacing_x) != 1) return 1;

    printf("Enter spacing Y in base units (e.g., 2): ");
    if (scanf("%d", &base_spacing_y) != 1) return 1;

    printf("Enter number of perturbed patterns: ");
    if (scanf("%d", &num_random_patterns) != 1) return 1;

    srand((unsigned int)time(NULL));

    // --- Porosity Calculation ---
    // Physical obstacle unit size = 1x1
    double obstacle_area = 1.0 * 1.0; 
    double unit_cell_area = (1.0 + base_spacing_x) * (1.0 + base_spacing_y);
    double porosity = 1.0 - (obstacle_area / unit_cell_area);
    double porosity_pct = porosity * 100.0;

    // Convert base dimensions to scaled sub-pixel dimensions
    int N_grid = N_base * scale_factor;
    int obj_size = 1 * scale_factor;              // Base 1x1 obstacle -> (scale_factor) x (scale_factor)
    int spacing_x = base_spacing_x * scale_factor;
    int spacing_y = base_spacing_y * scale_factor;
    int PAD = 5 * scale_factor;
    int shift_step = 1;                          // Moves by 1 sub-pixel unit

    int step_x = obj_size + spacing_x;
    int step_y = obj_size + spacing_y;
    int N_padded = N_grid + (2 * PAD);

    // Dynamic directory string incorporating porosity %
    char dir_name[128];
    sprintf(dir_name, "n%d_por_%.1f_sx%d_sy%d_scale%dx", 
            N_base, porosity_pct, base_spacing_x, base_spacing_y, scale_factor);
    MAKE_DIR(dir_name);

    // Allocate oversized buffer to accommodate extra edge obstacles
    Obstacle *obstacles = malloc((N_padded * N_padded) * sizeof(Obstacle));
    int *final_grid = (int *)calloc(N_grid * N_grid, sizeof(int));

    // ==========================================
    // STEP 3: Pattern 1 - Inline Grid Creation
    // (Shifted down by half an obstacle size to cover top/bottom channels)
    // ==========================================
    int center_shift = obj_size / 2;
    int obs_count = 0;

    for (int r = -center_shift; r <= N_padded; r += step_y) {
        for (int c = -center_shift; c <= N_padded; c += step_x) {
            obstacles[obs_count].r1 = r - PAD;
            obstacles[obs_count].r2 = r + obj_size - 1 - PAD;
            obstacles[obs_count].c1 = c - PAD;
            obstacles[obs_count].c2 = c + obj_size - 1 - PAD;
            obs_count++;
        }
    }

    render_to_grid(obstacles, obs_count, final_grid, N_grid);
    save_pattern(final_grid, N_grid, dir_name, "pattern_1");

    // ==========================================
    // STEP 4: Pattern 2 - Perfect Staggered Grid
    // ==========================================
    int vertical_offset = step_y / 2;
    Obstacle *staggered_obstacles = malloc(obs_count * sizeof(Obstacle));
    memcpy(staggered_obstacles, obstacles, obs_count * sizeof(Obstacle));

    for (int k = 0; k < obs_count; k++) {
        int col_idx = (staggered_obstacles[k].c1 + PAD + center_shift) / step_x;
        if (col_idx % 2 == 1) {
            staggered_obstacles[k].r1 += vertical_offset;
            staggered_obstacles[k].r2 += vertical_offset;
        }
    }

    render_to_grid(staggered_obstacles, obs_count, final_grid, N_grid);
    save_pattern(final_grid, N_grid, dir_name, "pattern_2");

    // ==========================================
    // STEP 5: Random Perturbation Patterns
    // ==========================================
    Obstacle *perturbed = malloc(obs_count * sizeof(Obstacle));

    for (int p = 0; p < num_random_patterns; p++) {
        memcpy(perturbed, staggered_obstacles, obs_count * sizeof(Obstacle));

        for (int k = 0; k < obs_count; k++) {
            int move = rand() % 5; // 0: Stay, 1: Up, 2: Down, 3: Left, 4: Right
            int dr = (move == 1) ? -shift_step : (move == 2) ? shift_step : 0;
            int dc = (move == 3) ? -shift_step : (move == 4) ? shift_step : 0;

            perturbed[k].r1 += dr;
            perturbed[k].r2 += dr;
            perturbed[k].c1 += dc;
            perturbed[k].c2 += dc;
        }

        render_to_grid(perturbed, obs_count, final_grid, N_grid);

        char pat_name[64];
        sprintf(pat_name, "pattern_%d", p + 3);
        save_pattern(final_grid, N_grid, dir_name, pat_name);
    }

    // Free memory
    free(obstacles);
    free(staggered_obstacles);
    free(perturbed);
    free(final_grid);

    printf("\nDone! Calculated Porosity: %.2f%%\nGenerated output in directory: '%s'\n", 
           porosity_pct, dir_name);

    return 0;
}
