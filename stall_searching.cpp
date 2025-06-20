#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iostream>
#include <string>
#include <unordered_set>
#include <omp.h>

using namespace std;

/*
*  Create all possible column configurations and filter based on fast conditions before storage.
*  A "config" is a vector of 'n' unsigned integers representing how many vertices in each column is unfilled.
*  Filters out configurations with incorrect number of filled vertices, fail basic column sum tests, or are reverses of existing configurations.
*/
vector<vector<unsigned>> make_configs(const vector<unsigned> dim, const int filled) {
    const int m = dim[0];
    const int n = dim[1];
    vector<vector<unsigned>> out = vector<vector<unsigned>>();
    // Each column can have up to m filled vertices, leading to m + 1 ^ n possible options for the entire graph
    unsigned long long int total_configs = pow(m+1, n);
    for (unsigned long long int i = 0; i < total_configs; i++) { //TODO parallelize this, but theres a race condition for checking reverses
        int config_sum = 0;
        vector<unsigned> config = vector<unsigned>(n);
        // Check whether the end columns are filled
        if (i % (m+1) == 0 || (int)floor(i / pow(m+1, n - 1)) % (m+1) == 0) {
            continue;
        }


        bool flag = false;
        // The last column with 0 unfilled
        int last_zero = 0;
        int one_counter = 0;
        //create the config
        for (int j = 0; j < n; j++) {
            // Calculate how many unfilled in the jth column
            config[j] = (int)floor(i / pow(m+1, j)) % (m+1);
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
        vector<unsigned> rev(n); //maybe do this at the end to remove race condition
        reverse_copy(config.begin(), config.end(), rev.begin());
        if (out.size() == 0 || find(out.begin(), out.end(), rev) == out.end()) {
            out.push_back(config);
        }
    }
    return out;

}

/* 
* Returns the first neighborhood indices of a particular index as a vector, with the first element being the original index. 
* All neighborhoods must be length 5 for the stall checking algorithm.
* Indices are COLUMN MAJOR!!!!
*/
vector<unsigned> find_nb(const unsigned index, const vector<unsigned> dim) {
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
vector<unsigned> find_hats(const unsigned index, const vector<unsigned> dim) {
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
vector<unsigned> find_hats_grid(const vector<unsigned> dim) {
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
* Checks whether a collection of zeroes stalls the grid according to the hats for the mxn grid.
* Any optimization to this function has large consequences since this is run many many times.
* Coming back, we can probably store the zeroes as a list and check the respective index since the hat vector is stored in order.
* Right now, we iterate over all of the indices and check whether its a 0 or a 1 iteratively.
*/
bool stall_check(unordered_set<unsigned> zeroes, vector<unsigned> hats) {
    vector<unsigned>::iterator i = hats.begin();
    while (i != hats.end()) {
        // if there is a 1 at the index, move onto the next index
        if (!zeroes.contains(*i)) { advance(i, 5); continue; }

        // *i == 0, each iteration is ONE hat
        vector<unsigned>::iterator next_check = next(i, 5);
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

int fact(int n)
{
    if (n == 0)
        return 1;
    int res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}

int nCk(int n, int k)
{
    return fact(n) / (fact(k) * fact(n - k));
}


//from rosetta code. idk how this works ngl.
vector<vector<unsigned>> find_combn(const int N, const int K, const int c, const vector<unsigned> dim) {
    vector<vector<unsigned>> combs(nCk(N, K));
    string bitmask(K, 1);
    bitmask.resize(N, 0);
    int j = 0;
    auto it = combs.begin();
    do {
        vector<unsigned> comb = vector<unsigned>(K);
        int k = 0;
        for (unsigned i = 0; i < N; ++i) {
            if (bitmask[i]) {
                comb[k] = i + 1 + (c - 1) * dim[0];
                ++k;
            }

        }
        *it = comb;
        ++it;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return combs;
}

/*
* Finds all stalled grids of a specific column configuration.
*/
vector<unordered_set<unsigned>> check_config_rev(const vector<unsigned> config,
    const vector<unsigned> dim,
    const vector<unsigned> hats,
    const unsigned last_zero) {
    const int m = dim[0];
    const int n = dim[1];

    vector<vector<vector<unsigned>>> options(n);
    vector<unsigned> lengths(n);
    vector<unsigned> factors(n);
    int out_length = accumulate(config.begin(), config.end(), 0);

    vector<unordered_set<unsigned>> stallout;
    ;
    //create indices
    for (int c = 1; c <= n; c++) {
        // All ways to choose which indices in the column are zeroes
        vector<vector<unsigned>> colIndices = find_combn(m, config[c - 1], c + (last_zero + 1), dim);

        options[c - 1] = colIndices;
        lengths[c - 1] = colIndices.size();
        if (c == 1) {
            factors[c - 1] = 1;
        }
        else {
            factors[c - 1] = lengths[c - 2] * factors[c - 2];
        }
    }

    int total_length = factors[n - 1] * lengths[n - 1];

    for (int i = 0; i < total_length; i++) {
        unordered_set<unsigned> zeroes(out_length);
        int zero_index = 0;
        for (int k = 0; k < n; k++) {
            vector<unsigned> option = options[k][i / factors[k] % lengths[k]];
            zeroes.insert(option.begin(), option.end());
        }

        if (stall_check(zeroes, hats)) {
            stallout.push_back(zeroes);
        }
    }


    return stallout;
}

vector<unsigned> create_reverse(const vector<unsigned> combn, const vector<unsigned> dim) {
    int n = combn.size();
    vector<unsigned> out(n);
    for (int i = 0; i < n; i++) {
        out[i] = (dim[0] + 1) - combn[i];
    }
    sort(out.begin(), out.end());
    return out;
}
vector<unordered_set<unsigned>> check_config(const vector<unsigned> config, const vector<unsigned> dim, const vector<unsigned> hats) {
    const int m = dim[0];
    const int n = dim[1];

    vector<vector<vector<unsigned>>> options(n);
    vector<int> lengths(n);
    vector<int> factors(n);
    int out_length = accumulate(config.begin(), config.end(), 0);

    vector<unordered_set<unsigned>> stallout;
    //create indices
    for (int c = 1; c <= n; c++) {
        vector<vector<unsigned>> colIndices = find_combn(m, config[c - 1], c, dim);

        int n_ = colIndices.size();
        if (c == 1) {

            // is the reverse in here already?
            vector<vector<unsigned>> noRevIndices = vector<vector<unsigned>>();
            noRevIndices.push_back(colIndices[0]);
            for (int r = 1; r < n_; r++) {
                if (find(noRevIndices.begin(), noRevIndices.end(), create_reverse(colIndices[r], dim)) == noRevIndices.end()) {
                    noRevIndices.push_back(colIndices[r]);
                }
            }

            options[c - 1] = noRevIndices;
            lengths[c - 1] = noRevIndices.size();
            factors[c - 1] = 1;

        }
        else {
            options[c - 1] = colIndices;
            lengths[c - 1] = colIndices.size();
            factors[c - 1] = lengths[c - 2] * factors[c - 2];
        }

    }

    int total_length = factors[dim[1] - 1] * lengths[dim[1] - 1];

    //iterate over all possible column configurations, check if one stalls
    for (int i = 0; i < total_length; i++) {
        unordered_set<unsigned> zeroes(out_length); //this is by far the slowest part of the algorithm
        int zero_index = 0;
        for (int k = 0; k < n; k++) {
            vector<unsigned> option = options[k][i / factors[k] % lengths[k]];
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
vector<unordered_set<unsigned>> check_config_spl(const vector<unsigned> config, const vector<unsigned> dim, const vector<unsigned> hats) {
    if (find(config.begin(), config.end(), 0) == config.end()) {
        return check_config(config, dim, hats);
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

    for (int j = 0; j <= dim[1]; j++) {
        if (j == dim[1] || config[j] == 0) {
            vector<unsigned> sub_dim = { dim[0], subcols };
            vector<unsigned> sub_hats = find_hats_grid(sub_dim);
            vector<unsigned> sub_config = vector<unsigned>(subcols);
            copy_n(config.begin() + (last_zero + 1), subcols, sub_config.begin());
            vector<unordered_set<unsigned>> sub_stall = check_config_rev(sub_config, sub_dim, sub_hats, last_zero);

            if (sub_stall.size() == 0) {
                break;
            }

            sub_stalls[which_sub] = sub_stall;
            lengths[which_sub] = sub_stall.size();

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
        int zero_index = 0;
        zeroes.reserve(out_length);
        for (int k = 0; k < num_subs; k++) {
            unordered_set<unsigned> option = sub_stalls[k][i / factors[k] % lengths[k]];
            zeroes.insert(option.begin(), option.end());
        }

        if (stall_check(zeroes, hats)) {
            stallout.push_back(zeroes);
        }
    }
    return stallout;


}

vector<unordered_set<unsigned>> stall_search(const vector<unsigned> dim, vector<vector<unsigned>> configs) {
    vector<unsigned> hats = find_hats_grid(dim);
    vector<unordered_set<unsigned>> out;
    printf("Number of configs: %u \nProgress: ", configs.size());
    for (auto it = configs.begin(); it != configs.end(); it++) {
        vector<unordered_set<unsigned>> stalls = check_config_spl(*it, dim, hats);
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

vector<unordered_set<unsigned>> stall_search_par(const vector<unsigned> dim, const vector<vector<unsigned>> configs) {
    vector<unsigned> hats = find_hats_grid(dim);
    vector<unordered_set<unsigned>> out;
    printf("Number of configs: %u \nProgress: ", configs.size());
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < configs.size(); i++) {
        vector<unordered_set<unsigned>> stalls = check_config_spl(configs[i], dim, hats);
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



void print_grids(vector<unordered_set<unsigned>> zeroes, const vector<unsigned> dim) {
    int num_grids = zeroes.size();
    vector<unsigned> grid(dim[0] * dim[1], 1);
    unordered_set<unsigned> option;
    for (int k = 0; k < num_grids; k++) {
        option = zeroes[k];
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
