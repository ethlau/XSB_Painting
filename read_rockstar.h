#define BINARY_HEADER_SIZE 256
#define VERSION_MAX_SIZE 12
#define TOTAL_CHUNK 216

struct binary_output_header {
    unsigned long long magic;
    long long snap, chunk;
    float scale, Om, Ol, h0;
    float bounds[6];
    long long num_halos, num_particles;
    float box_size, particle_mass;
    long long particle_type;
    int format_revision;
    char rockstar_version[VERSION_MAX_SIZE];
    char unused[BINARY_HEADER_SIZE - (sizeof(char)*VERSION_MAX_SIZE) - (sizeof(float)*12) - sizeof(int) - (sizeof(long long)*6)];
};

struct halo {
    long long id;
    float pos[6], corevel[3], bulkvel[3];
    float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
          J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
          b_to_a2, c_to_a2, A2[3],
          bullock_spin, kin_to_pot, m_pe_b, m_pe_d, halfmass_radius;
    long long num_p, num_child_particles, p_start, desc, flags, n_core;
    float min_pos_err, min_vel_err, min_bulkvel_err;
};

void load_rockstar(const char *root, int snap, struct binary_output_header *bheader, struct halo **halos, long long **part_ids, float **part_pos, long long *TotHalo, long long *TotIds);
void load_header_rockstar(const char *root, int snap, int chunk, struct binary_output_header *bheader);
void load_halos_rockstar(const char *root, int snap, int chunk, struct binary_output_header *bheader, struct halo **halos, long long **part_ids, float **part_pos, long long nhalo, long long nids);
void load_parent(const char *root, int snap, struct halo *halos, long long tothalos);
