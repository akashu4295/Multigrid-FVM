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
    int r1, r2, c1, c2;
} Obstacle;

int main() {
    int N_grid;
    double ratio;
    int num_patterns;
    double pct1, pct2, pct3;

    printf("Enter N_grid (e.g., 50 for a 50x50 grid): ");
    if (scanf("%d", &N_grid) != 1) return 1;

    if (N_grid < 7) {
        printf("N_grid must be at least 7 to accommodate 3x3 particles safely within bounds.\n");
        return 1;
    }

    printf("Enter ratio of grids to be filled (e.g., 0.3 for 30%%): ");
    if (scanf("%lf", &ratio) != 1) return 1;

    printf("Enter number of patterns to create: ");
    if (scanf("%d", &num_patterns) != 1) return 1;

    printf("Enter percentage of total particles that are size 1 (1x1) (e.g., 50): ");
    if (scanf("%lf", &pct1) != 1) return 1;
    
    printf("Enter percentage of total particles that are size 2 (2x2) (e.g., 30): ");
    if (scanf("%lf", &pct2) != 1) return 1;
    
    printf("Enter percentage of total particles that are size 3 (3x3) (e.g., 20): ");
    if (scanf("%lf", &pct3) != 1) return 1;

    if (pct1 + pct2 + pct3 > 100.1 || pct1 + pct2 + pct3 < 99.9) {
        printf("Error: Percentages must sum to 100.\n");
        return 1;
    }

    srand((unsigned int)time(NULL)); // seed point based on time

    char dir_name[100];
    sprintf(dir_name, "n_%d/p_%.2f", N_grid, ratio);
    MAKE_DIR(dir_name);
    printf("Saving outputs to folder: %s\n", dir_name);

    int *grid = (int *)malloc(N_grid * N_grid * sizeof(int));
    if (grid == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    double p1 = pct1 / 100.0;
    double p2 = pct2 / 100.0;
    double p3 = pct3 / 100.0;

    int valid_rows = N_grid; // Rows 1 to N_grid
    int valid_cols = N_grid; // Cols 1 to N_grid
    
    double target_area = valid_rows * valid_cols * ratio;
    double avg_area_per_particle = (p1 * 1.0) + (p2 * 4.0) + (p3 * 9.0);
    int total_particles = (int)(target_area / avg_area_per_particle + 0.5);

    int target_p1 = (int)(total_particles * p1 + 0.5);
    int target_p2 = (int)(total_particles * p2 + 0.5);
    int target_p3 = (int)(total_particles * p3 + 0.5);

    printf("\nTargeting ~%d total particles to maintain ratio.\n", total_particles);
    printf("Size 1: %d | Size 2: %d | Size 3: %d\n\n", target_p1, target_p2, target_p3);

    for (int p = 0; p < num_patterns; p++) {
        memset(grid, 0, N_grid * N_grid * sizeof(int));
        
        Obstacle *obstacles = malloc((target_p1 + target_p2 + target_p3) * sizeof(Obstacle));
        int obs_count = 0;

        int sizes[] = {3, 2, 1};
        int counts[] = {target_p3, target_p2, target_p1};

        for (int s = 0; s < 3; s++) {
            int size = sizes[s];
            int needed = counts[s];

            for (int k = 0; k < needed; k++) {
                int placed = 0;
                int max_retries = 10000; // Prevent infinite loop if density is impossibly high

                for (int ret = 0; ret < max_retries; ret++) {
                    int r = rand() % (valid_rows - size + 1);
                    int c = rand() % (valid_cols - size + 1);

                    // Check for overlap
                    int overlap = 0;
                    for (int i = 0; i < size; i++) {
                        for (int j = 0; j < size; j++) {
                            if (grid[(r + i) * N_grid + (c + j)] == 1) {
                                overlap = 1;
                                break;
                            }
                        }
                        if (overlap) break;
                    }

                    // Place the particle if space is free
                    if (!overlap) {
                        for (int i = 0; i < size; i++) {
                            for (int j = 0; j < size; j++) {
                                grid[(r + i) * N_grid + (c + j)] = 1;
                            }
                        }
                        
                        obstacles[obs_count].r1 = r+1;
                        obstacles[obs_count].r2 = r + size;
                        obstacles[obs_count].c1 = c+1;
                        obstacles[obs_count].c2 = c + size;
                        obs_count++;
                        
                        placed = 1;
                        break;
                    }
                }
                if (!placed) {
                    printf("Warning in pattern %d: Could not fit a size %dx%d particle. Density may be too high.\n", p+1, size, size);
                }
            }
        }

        char txt_filename[256];
        sprintf(txt_filename, "%s/pattern_%d.txt", dir_name, p + 1);
        FILE *txt_file = fopen(txt_filename, "w");
        
        if (txt_file) {
            fprintf(txt_file, "%d\n", obs_count);
            for (int i = 0; i < obs_count; i++) {
                // Formatting: row_start row_end col_start col_end
                fprintf(txt_file, "%d %d %d %d\n", obstacles[i].r1, obstacles[i].r2, obstacles[i].c1, obstacles[i].c2);
            }
            fclose(txt_file);
        } else {
            printf("Error creating file: %s\n", txt_filename);
        }

        char img_filename[256];
        sprintf(img_filename, "%s/pattern_%d.pgm", dir_name, p + 1);
        FILE *img_file = fopen(img_filename, "w");
        
        if (img_file) {
            fprintf(img_file, "P2\n%d %d\n255\n", N_grid, N_grid);
            
            for (int i = 0; i < N_grid; i++) {
                for (int j = 0; j < N_grid; j++) {
                    int color = (grid[i * N_grid + j] == 1) ? 0 : 255;
                    fprintf(img_file, "%d ", color);
                }
                fprintf(img_file, "\n");
            }
            fclose(img_file);
        } else {
            printf("Error creating file: %s\n", img_filename);
        }

        free(obstacles);
    }

    free(grid);
    printf("Successfully generated %d patterns.\n", num_patterns);

    return 0;
}
