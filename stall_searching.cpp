#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <string>
#include <unordered_set>
#include <set>
#include <atomic>
#include <omp.h>

using namespace std;

/*
*  Create all possible column configurations and filter based on fast conditions before storage.
*  A "config" is a vector of 'n' unsigned integers representing how many vertices in each column is unfilled.
*  Filters out configurations with incorrect number of filled vertices, fail basic column sum tests, or are reverses of existing configurations.
*/
vector<vector<unsigned>> make_configs(const vector<unsigned>& dim, const int filled) {
    const int m = dim[0];
    const int n = dim[1];
    vector<vector<unsigned>> out;

    // Precompute powers of (m+1) once so the mixed-radix digit decomposition
    // below is plain integer division instead of pow()/floor() every
    // iteration (floating point pow() also silently loses precision once
    // total_configs exceeds 2^53).
    vector<unsigned long long> pow_table(n + 1);
    pow_table[0] = 1;
    for (int j = 1; j <= n; j++) {
        pow_table[j] = pow_table[j - 1] * (unsigned long long)(m + 1);
    }
    // Each column can have up to m filled vertices, leading to m + 1 ^ n possible options for the entire graph
    unsigned long long total_configs = pow_table[n];

    // Mirrors the contents of `out` so the "is my reverse already accepted"
    // check below is O(log |out|) instead of the O(|out|) linear scan it
    // replaces -- the original made this function O(|out|^2) overall.
    set<vector<unsigned>> seen;

    for (unsigned long long i = 0; i < total_configs; i++) {
        //TODO parallelize this -- the reverse-dedup below is now a hash/tree lookup (seen), so the former
        // concern about scanning `out` is gone, but accept/reject still needs to be serialized against
        // `seen`/`out` (e.g. a critical section, or per-thread buffers merged and re-deduped at the end).
        int config_sum = 0;
        vector<unsigned> config(n);
        // Check whether the end columns are filled
        if (i % (m + 1) == 0 || (i / pow_table[n - 1]) % (m + 1) == 0) {
            continue;
        }

        bool flag = false;
        // The last column with 0 unfilled
        int last_zero = 0;
        int one_counter = 0;
        //create the config
        for (int j = 0; j < n; j++) {
            // Calculate how many unfilled in the jth column
            config[j] = (unsigned)((i / pow_table[j]) % (m + 1));
            config_sum += config[j];

            // If we have too many unfilled, continue
            if (config_sum > m * n - filled) {
                flag = true;
                break;
            }

            // If the second (or second to last) column is full, the first column cannot force as a path. Special case of column splitting.
            // If I add dynamic solutions, I can expland this to help solve general grid splitting problems
            if (j > 0) {
                if (config[1] == 0 && config[0] < ceil((float)(m + 1) / 2)) {
                    flag = true;
                    break;
                }
            }

            // If the last column had no unfilled vertices, this column must have the same number of unfilled vertices as the column opposite the fully filled column.
            // In general, it must also have the same rows unfilled as the column opposite, but we can't check that now.
            if (j > 1 && last_zero == j - 1 &&
                config[j - 2] != config[j]) {
                flag = true;
                break;
            }

            if (config[j] == 0) {
                // This condition is a bit odd. Consider a sequence of columns with only one vertex unfilled sandwiched between two filled columns.
                // If the number of columns is less than m, then that block does not stall.
                // This can be removed in a dynamic solution algorithm since this block would have more vertices than a maximum FZF set.
                if (one_counter < m) {
                    flag = true;
                    break;
                }
                one_counter = 0;

                if (j > 1) {
                    // 2 adjacent filled columns fill a grid xd
                    if (last_zero == j - 1) {
                        flag = true;
                        break;
                    }

                    // More grid splitting xd.
                    if (j > 2) {
                        if (config[j - 1] < ceil((float)(m + 1) / 2) &&
                            last_zero == j - 2) {
                            flag = true;
                            break;
                        }
                    }
                }

                last_zero = j;
            }
            else if (config[j] == 1) {
                one_counter++;
            }
            else {
                one_counter += m;
            }

        }
        if (flag || one_counter < m) {
            continue;
        }

        if (config[n - 2] == 0 && config[n - 1] < ceil((float)(m + 1) / 2)) {
            continue;
        }

        //check sum of config
        if (config_sum < m * n - filled) {
            continue;
        }

        // Check if a reversed copy of the column configuration is already stalling.
        vector<unsigned> rev(n);
        reverse_copy(config.begin(), config.end(), rev.begin());
        if (seen.find(rev) == seen.end()) {
            out.push_back(config);
            seen.insert(config);
        }
    }
    return out;

}

/*
* Returns the first neighborhood indices of a particular index as a vector, with the first element being the original index.
* All neighborhoods must be length 5 for the stall checking algorithm.
* Indices are COLUMN MAJOR!!!!
*/
vector<unsigned> find_nb(const unsigned index, const vector<unsigned>& dim) {
    const int m = dim[0];
    const int n = dim[1];
    vector<unsigned> nb;

    // Path special case.
    if (n == 1) {
        if (index == 1) {
            nb = { index, index + 1, 0, 0, 0 };
        }
        else if (index == m) {
            nb = { index, index - 1, 0, 0, 0 };
        }
        else {
            nb = { index, index + 1, index - 1, 0, 0 };
        }

        return nb;
    }
    // Top row
    if (index % m == 1) {
        if (index <= m) {
            nb = { index, index + 1, index + m, 0, 0 };
        }
        else if (index > (n * m - m)) {
            nb = { index, index + 1, index - m, 0, 0 };
        }
        else {
            nb = { index, index + 1, index - m, index + m, 0 };
        }
    }
    else if (index % m == 0) { // Bottom row
        if (index <= m) {
            nb = { index, index - 1, index + m, 0, 0 };
        }
        else if (index > (n * m - m)) {
            nb = { index, index - 1, index - m, 0, 0 };
        }
        else {
            nb = { index, index - 1, index - m, index + m, 0 };
        }
    }
    else {
        if (index <= m) { // First column
            nb = { index, index - 1, index + 1, index + m, 0 };
        }
        else if (index > (n * m - m)) { // Last column
            nb = { index, index - 1, index + 1, index - m, 0 };
        }
        else {
            nb = { index, index - 1, index + 1, index - m, index + m };
        }
    }

    return nb;
}

/*
* Returns the hat (second neighborhood) indices for a given index.
* Because we assume there are more filled vertices than unfilled,
* checking the hats for the unfilled vertices is WAY faster than checking the forcing condition for all filled vertices.
*/
vector<unsigned> find_hats(const unsigned index, const vector<unsigned>& dim) {
    vector<unsigned> ind_nb = find_nb(index, dim);
    unsigned num_neighbors = std::count_if(ind_nb.begin() + 1, ind_nb.end(), [](unsigned i) {return i > 0; });
    vector<unsigned> hats(5 * num_neighbors);
    unsigned hat_index = 0;
    for (unsigned i = 1; i <= num_neighbors; i++) {
        unsigned neighbor = ind_nb[i];
        vector<unsigned> nb_hat = find_nb(neighbor, dim);
        iter_swap(nb_hat.begin(), find(nb_hat.begin(), nb_hat.end(), index));
        std::copy(nb_hat.begin(), nb_hat.end(), hats.begin() + hat_index);
        hat_index += 5;
    }

    return hats;
}

/*
* Find the second neighborhoods for the entire grid.
* Precalculating all of these since it gets checked millions of times.
*/
vector<unsigned> find_hats_grid(const vector<unsigned>& dim) {
    const unsigned m = dim[0];
    const unsigned n = dim[1];

    vector<unsigned> hats;
    if (n == 1) {
        hats = vector<unsigned>(5 * (2 * 1 + 2 * (m - 2)));

    }
    else {
        hats = vector<unsigned>(5 * (2 * 4 + 6 * (m - 2 + n - 2) + 4 * ((m - 2) * (n - 2))));
    }
    unsigned hat_it = 0;

    for (unsigned i = 0; i < m * n; i++) {
        vector<unsigned> ind_nb = find_nb(i + 1, dim);
        unsigned num_neighbors = std::count_if(ind_nb.begin(), ind_nb.end(), [](unsigned i) {return i > 0; }) - 1;
        vector<unsigned> index_hats = find_hats(i + 1, dim);
        std::copy(index_hats.begin(), index_hats.end(), hats.begin() + hat_it);
        hat_it += 5 * num_neighbors;
    }

    return hats;
}

/*
* Precomputes find_hats_grid for every possible sub-grid column count (1..n) once per run.
* check_config_spl needs the hats for shape {m, subcols} for whatever subcols a column split
* produces, and the same shapes recur across many configs -- build the table once instead of
* recomputing it from scratch for every config that happens to split the same way.
*/
vector<vector<unsigned>> precompute_subcol_hats(const unsigned m, const unsigned n) {
    vector<vector<unsigned>> table(n + 1);
    for (unsigned subcols = 1; subcols <= n; subcols++) {
        table[subcols] = find_hats_grid({ m, subcols });
    }
    return table;
}

/*
* Checks whether a collection of zeroes stalls the grid according to the hats for the mxn grid.
* Any optimization to this function has large consequences since this is run many many times.
* Coming back, we can probably store the zeroes as a list and check the respective index since the hat vector is stored in order.
* Right now, we iterate over all of the indices and check whether its a 0 or a 1 iteratively.
*/
bool stall_check(const unordered_set<unsigned>& zeroes, const vector<unsigned>& hats) {
    vector<unsigned>::const_iterator i = hats.begin();
    while (i != hats.end()) {
        // if there is a 1 at the index, move onto the next index
        if (!zeroes.contains(*i)) { advance(i, 5); continue; }

        // *i == 0, each iteration is ONE hat
        vector<unsigned>::const_iterator next_check = next(i, 5);
        bool flag = false;
        ++i;

        //check all nodes in the hat
        while (*i != 0 && i != next_check) {
            //if there is a zero in the hat, we are happy and move onto the next hat
            if (zeroes.contains(*i)) {
                flag = true;
                break;
            }
            ++i;
        }
        //the zero we found gets filled in the next step! the set doesn't stall :(
        if (!flag) {
            return false;
        }

        i = next_check;
    }
    return true;
}

// Overflow-safe iterative binomial coefficient (the previous fact()-based version overflowed
// 32-bit int around n=13, silently corrupting anything sized from it).
unsigned long long nCk(int n, int k)
{
    if (k < 0 || k > n) return 0;
    k = min(k, n - k);
    unsigned long long result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (unsigned long long)(n - i) / (unsigned long long)(i + 1);
    }
    return result;
}

/*
* Precomputes, for every possible per-column unfilled count k (0..m), the set of k-subsets of a
* single column's m rows -- expressed as row indices (1..m) with no column offset applied.
* Every config/column that needs k unfilled rows reuses this table via a cheap offset add
* (see offset_combo) instead of re-deriving the same combinatorics with prev_permutation from
* scratch for every column of every config.
*/
vector<vector<vector<unsigned>>> precompute_base_combos(const unsigned m) {
    vector<vector<vector<unsigned>>> table(m + 1);
    for (unsigned k = 0; k <= m; k++) {
        vector<vector<unsigned>> combos((size_t)nCk((int)m, (int)k));
        string bitmask(k, 1);
        bitmask.resize(m, 0);
        auto it = combos.begin();
        do {
            vector<unsigned> comb(k);
            unsigned c = 0;
            for (unsigned bi = 0; bi < m; ++bi) {
                if (bitmask[bi]) {
                    comb[c] = bi + 1;
                    ++c;
                }
            }
            *it = comb;
            ++it;
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
        table[k] = move(combos);
    }
    return table;
}

// Applies a column offset (in whole columns, i.e. rows-per-column units) to a cached base combo.
inline vector<unsigned> offset_combo(const vector<unsigned>& base, const unsigned col_offset_rows) {
    vector<unsigned> out(base.size());
    for (size_t i = 0; i < base.size(); i++) {
        out[i] = base[i] + col_offset_rows;
    }
    return out;
}

vector<unsigned> create_reverse(const vector<unsigned>& combn, const vector<unsigned>& dim) {
    int n = combn.size();
    vector<unsigned> out(n);
    for (int i = 0; i < n; i++) {
        out[i] = (dim[0] + 1) - combn[i];
    }
    sort(out.begin(), out.end());
    return out;
}

/*
* Finds all stalled grids of a specific column configuration.
*/
vector<unordered_set<unsigned>> check_config_rev(const vector<unsigned>& config,
    const vector<unsigned>& dim,
    const vector<unsigned>& hats,
    const unsigned last_zero,
    const vector<vector<vector<unsigned>>>& base_combos) {
    const int m = dim[0];
    const int n = dim[1];

    vector<vector<vector<unsigned>>> options(n);
    vector<unsigned> lengths(n);
    vector<unsigned> factors(n);
    int out_length = accumulate(config.begin(), config.end(), 0);

    vector<unordered_set<unsigned>> stallout;
    //create indices
    for (int c = 1; c <= n; c++) {
        // All ways to choose which indices in the column are zeroes
        const vector<vector<unsigned>>& base = base_combos[config[c - 1]];
        unsigned col_offset_rows = (unsigned)(c + (int)last_zero) * (unsigned)m;

        vector<vector<unsigned>> colIndices;
        colIndices.reserve(base.size());
        for (const auto& relCombo : base) {
            colIndices.push_back(offset_combo(relCombo, col_offset_rows));
        }

        options[c - 1] = move(colIndices);
        lengths[c - 1] = options[c - 1].size();
        if (c == 1) {
            factors[c - 1] = 1;
        }
        else {
            factors[c - 1] = lengths[c - 2] * factors[c - 2];
        }
    }

    int total_length = factors[n - 1] * lengths[n - 1];

    unordered_set<unsigned> zeroes;
    zeroes.reserve(out_length);
    for (int i = 0; i < total_length; i++) {
        zeroes.clear();
        for (int k = 0; k < n; k++) {
            const vector<unsigned>& option = options[k][i / factors[k] % lengths[k]];
            zeroes.insert(option.begin(), option.end());
        }

        if (stall_check(zeroes, hats)) {
            stallout.push_back(zeroes);
        }
    }


    return stallout;
}

vector<unordered_set<unsigned>> check_config(const vector<unsigned>& config, const vector<unsigned>& dim,
    const vector<unsigned>& hats, const vector<vector<vector<unsigned>>>& base_combos) {
    const int m = dim[0];
    const int n = dim[1];

    vector<vector<vector<unsigned>>> options(n);
    vector<int> lengths(n);
    vector<int> factors(n);
    int out_length = accumulate(config.begin(), config.end(), 0);

    vector<unordered_set<unsigned>> stallout;
    //create indices
    for (int c = 1; c <= n; c++) {
        const vector<vector<unsigned>>& base = base_combos[config[c - 1]];

        if (c == 1) {
            // is the reverse in here already?
            vector<vector<unsigned>> noRevIndices;
            set<vector<unsigned>> accepted; // mirrors noRevIndices for an O(log n) membership check
            for (size_t r = 0; r < base.size(); r++) {
                vector<unsigned> combo = offset_combo(base[r], 0);
                if (r == 0) {
                    noRevIndices.push_back(combo);
                    accepted.insert(combo);
                    continue;
                }
                vector<unsigned> rev = create_reverse(combo, dim);
                if (accepted.find(rev) == accepted.end()) {
                    accepted.insert(combo);
                    noRevIndices.push_back(move(combo));
                }
            }

            options[c - 1] = move(noRevIndices);
            lengths[c - 1] = options[c - 1].size();
            factors[c - 1] = 1;

        }
        else {
            vector<vector<unsigned>> colIndices;
            colIndices.reserve(base.size());
            unsigned col_offset_rows = (unsigned)(c - 1) * (unsigned)m;
            for (const auto& relCombo : base) {
                colIndices.push_back(offset_combo(relCombo, col_offset_rows));
            }
            options[c - 1] = move(colIndices);
            lengths[c - 1] = options[c - 1].size();
            factors[c - 1] = lengths[c - 2] * factors[c - 2];
        }

    }

    int total_length = factors[dim[1] - 1] * lengths[dim[1] - 1];

    unordered_set<unsigned> zeroes;
    zeroes.reserve(out_length);
    //iterate over all possible column configurations, check if one stalls
    for (int i = 0; i < total_length; i++) {
        zeroes.clear();
        for (int k = 0; k < n; k++) {
            const vector<unsigned>& option = options[k][i / factors[k] % lengths[k]];
            zeroes.insert(option.begin(), option.end());
        }

        if (stall_check(zeroes, hats)) {
            stallout.push_back(zeroes);
        }
    }


    return stallout;
}

// this handles the case where the grid can be split into multiple smaller grids.
// theres a dynamic programming solution for this that's way more elegant
// this function becomes astronomically slow for small m - large n scenarios but those can be trivially solved by hand
vector<unordered_set<unsigned>> check_config_spl(const vector<unsigned>& config, const vector<unsigned>& dim,
    const vector<unsigned>& hats,
    const vector<vector<unsigned>>& subcol_hats,
    const vector<vector<vector<unsigned>>>& base_combos) {
    if (find(config.begin(), config.end(), 0) == config.end()) {
        return check_config(config, dim, hats, base_combos);
    }

    unsigned subcols = 0; //preallocating a lot of values, for minor performance gains
    int which_sub = 0;
    int last_zero = -1;
    int num_subs = count_if(config.begin(), config.end(), [](unsigned u) {return u == 0; }) + 1;

    int out_length = accumulate(config.begin(), config.end(), 0);
    vector<unsigned> lengths(num_subs);
    vector<unsigned> factors(num_subs);
    unordered_set<unsigned> zeroes(out_length);

     //some of these could probably be converted to sets, but there isnt really a point
    vector<unordered_set<unsigned>> stallout;
    vector<vector<unordered_set<unsigned>>> sub_stalls(num_subs);

    for (int j = 0; j <= (int)dim[1]; j++) {
        if (j == (int)dim[1] || config[j] == 0) {
            vector<unsigned> sub_dim = { dim[0], subcols };
            const vector<unsigned>& sub_hats = subcol_hats[subcols];
            vector<unsigned> sub_config = vector<unsigned>(subcols);
            copy_n(config.begin() + (last_zero + 1), subcols, sub_config.begin());
            vector<unordered_set<unsigned>> sub_stall = check_config_rev(sub_config, sub_dim, sub_hats, last_zero, base_combos);

            if (sub_stall.size() == 0) {
                break;
            }

            sub_stalls[which_sub] = move(sub_stall);
            lengths[which_sub] = sub_stalls[which_sub].size();

            if (which_sub == 0) {
                factors[which_sub] = 1;
            }
            else {
                factors[which_sub] = lengths[which_sub - 1] * factors[which_sub - 1];
            }

            last_zero = j;
            which_sub++;
            subcols = 0;
        }
        else {
            subcols++;
        }
    }

    int total_length = factors[num_subs - 1] * lengths[num_subs - 1];

    for (unsigned i = 0; i < total_length; i++) {
        zeroes.clear();
        zeroes.reserve(out_length);
        for (int k = 0; k < num_subs; k++) {
            const unordered_set<unsigned>& option = sub_stalls[k][i / factors[k] % lengths[k]];
            zeroes.insert(option.begin(), option.end());
        }

        if (stall_check(zeroes, hats)) {
            stallout.push_back(zeroes);
        }
    }
    return stallout;


}

// Rough relative cost of enumerating a config's zero placements (product of per-column
// candidate counts, in log-space to avoid overflow). Used only to order the work queue so
// the heaviest configs get dispatched first under dynamic scheduling.
double estimate_config_log_cost(const vector<unsigned>& config, const vector<vector<vector<unsigned>>>& base_combos) {
    double log_cost = 0.0;
    for (unsigned k : config) {
        log_cost += log((double)base_combos[k].size() + 1.0);
    }
    return log_cost;
}

vector<unordered_set<unsigned>> stall_search(const vector<unsigned>& dim, const vector<vector<unsigned>>& configs) {
    vector<unsigned> hats = find_hats_grid(dim);
    vector<vector<unsigned>> subcol_hats = precompute_subcol_hats(dim[0], dim[1]);
    vector<vector<vector<unsigned>>> base_combos = precompute_base_combos(dim[0]);
    vector<unordered_set<unsigned>> out;
    printf("Number of configs: %zu \nProgress: ", configs.size());
    for (auto it = configs.begin(); it != configs.end(); it++) {
        vector<unordered_set<unsigned>> stalls = check_config_spl(*it, dim, hats, subcol_hats, base_combos);
        printf(".");
        if (stalls.size() > 0) {
            for (auto it2 = stalls.begin(); it2 != stalls.end(); it2++) {
                out.push_back(*it2);
                printf("!");
            }
        }
    }
    return out;
}

vector<unordered_set<unsigned>> stall_search_par(const vector<unsigned>& dim, const vector<vector<unsigned>>& configs) {
    vector<unsigned> hats = find_hats_grid(dim);
    vector<vector<unsigned>> subcol_hats = precompute_subcol_hats(dim[0], dim[1]);
    vector<vector<vector<unsigned>>> base_combos = precompute_base_combos(dim[0]);
    printf("Number of configs: %zu \nProgress: ", configs.size());

    // Longest-job-first: sort configs by estimated cost before handing them to dynamic
    // scheduling, so the few very expensive configs don't end up serialized at the tail.
    vector<int> order(configs.size());
    iota(order.begin(), order.end(), 0);
    vector<double> cost(configs.size());
    for (size_t i = 0; i < configs.size(); i++) {
        cost[i] = estimate_config_log_cost(configs[i], base_combos);
    }
    sort(order.begin(), order.end(), [&](int a, int b) { return cost[a] > cost[b]; });

    // Each thread accumulates into its own slot and results are merged after the parallel
    // region -- the original called push_back on a single shared vector from every thread,
    // which is a data race (concurrent reallocation can corrupt the vector for all threads).
    int num_threads = omp_get_max_threads();
    vector<vector<unordered_set<unsigned>>> thread_out(num_threads);
    atomic<size_t> progress_counter(0);
    size_t total_configs = configs.size();

    #pragma omp parallel for schedule(dynamic)
    for (int oi = 0; oi < (int)order.size(); oi++) {
        int i = order[oi];
        int tid = omp_get_thread_num();
        vector<unordered_set<unsigned>> stalls = check_config_spl(configs[i], dim, hats, subcol_hats, base_combos);
        if (!stalls.empty()) {
            thread_out[tid].insert(thread_out[tid].end(),
                make_move_iterator(stalls.begin()), make_move_iterator(stalls.end()));
        }

        // Progress is reported periodically by whichever thread crosses a multiple of 100,
        // instead of once per config/stall -- printf takes an internal stdio lock, so calling
        // it on every single config serialized the threads against each other.
        size_t done = ++progress_counter;
        if (done % 100 == 0 || done == total_configs) {
            #pragma omp critical
            {
                printf("\rProgress: %zu/%zu", done, total_configs);
                fflush(stdout);
            }
        }
    }
    printf("\n");

    vector<unordered_set<unsigned>> out;
    size_t total_out = 0;
    for (auto& v : thread_out) total_out += v.size();
    out.reserve(total_out);
    for (auto& v : thread_out) {
        out.insert(out.end(), make_move_iterator(v.begin()), make_move_iterator(v.end()));
    }
    return out;
}



void print_grids(const vector<unordered_set<unsigned>>& zeroes, const vector<unsigned>& dim) {
    int num_grids = zeroes.size();
    vector<unsigned> grid(dim[0] * dim[1], 1);
    for (int k = 0; k < num_grids; k++) {
        const unordered_set<unsigned>& option = zeroes[k];
        grid.assign(dim[0] * dim[1], 1);
        printf("Zero Indices: ");
        for (auto i = option.begin(); i != option.end(); i++) {
            printf("%u ", *i);
            grid[*i - 1] = 0;
        }
        printf("\n");
        printf("Grid:\n");
        for (int j = 0; j < dim[1]; j++) {
            for (int i = 0; i < dim[0]; i++) {
                printf("%u ", grid[j * dim[0] + i]);
            }
            printf("\n");
        }
        printf("\n");
        printf("\n");
    }
}

int main()
{
    unsigned m;
    unsigned n;
    int filled;
    int t;
    cout << "Search for Stalled Grids \n";
    cout << "Rows:";
    cin >> m;
    cout << "Columns:";
    cin >> n;
    cout << "Filled Vertices:";
    cin >> filled;
    cout << "Threads:";
    cin >> t;
    printf("\n");
    const vector<unsigned> dim = { m,n };
    omp_set_num_threads(t);


    clock_t start = clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    vector<vector<unsigned>> cols = make_configs(dim, filled);


    vector<unordered_set<unsigned>> stalls = stall_search_par(dim, cols);

    clock_t end = clock();
    auto t_end = std::chrono::high_resolution_clock::now();

    printf("\n\n");
    print_grids(stalls, dim);

    cout << "Dimension: (" << m << "," << n << "), Filled: " << filled << "\n";
    cout
        << "Wall clock time passed: "
        << std::chrono::duration<double, std::milli>(t_end - t_start).count() << " ms; "
        << std::chrono::duration<double>(t_end - t_start).count() << " s \n\n";

    system("pause");
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
