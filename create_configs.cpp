#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

static vector<uvec> make_configs(const uvec dim, const int filled) {
	const int m = dim[1];
	const int n = dim[2];
	vector<uvec> out;
	unsigned long long int total_configs = pow(m, n);
	//uvec config(n);
	//for (unsigned long long int i = 0; i < total_configs; i++) {
	//	//fast to check conditions
	//	if (i % m == 0 || i % (int)pow(m,2) || (int)floor(i / pow(m, n-1)) % m == 0) {
	//		continue;
	//	}

	//	//create the config
	//	bool flag = false;
	//	int last_zero = 0;
	//	for (int j = 0; j < n; j++) {
	//		config[j] = (int)floor(i / pow(m, j)) % m;
	//		if (j > 0) {
	//			if (config[1] == 0 && config[0] < std::ceil((m + 1) / 2)) {
	//				flag = true;
	//				break;
	//			}
	//		}

	//		if (j > 1 && last_zero == j - 1 &&
	//			config[j-2] != config[j]) {
	//			flag = true;
	//			break;
	//		}

	//		if (config[j] == 0) {
	//			if (j > 1) {
	//				if (last_zero == j - 1) {
	//					flag = true;
	//					break;
	//				}

	//				if (j > 2) {
	//					if (config[j - 1] < std::ceil((m + 1) / 2) &&
	//						last_zero == j - 2) {
	//						flag = true;
	//						break;
	//					}
	//				}
	//			}

	//			last_zero = j;
	//		}
	//		
	//	}
	//	if (flag) {
	//		continue;
	//	}

	//	if (config[n - 1] == 0 && config[n - 2] < std::ceil((m + 1) / 2)) {
	//		continue;
	//	}

	//	//check sum of config
	//	if (accu(config) != m * n - filled) {
	//		continue;
	//	}

	//	//harder to check conditions
	//	if (any(config == 0)) {
	//		int one_counter = 0;
	//		for (int j = 0; j < n; j++) {
	//			if (config[j] == 0) {
	//				if (one_counter < m) {
	//					flag = true;
	//				}
	//			}
	//			else if (config[j] == 1) {
	//				one_counter++;
	//			}
	//			else {
	//				one_counter += m;
	//			}
	//		}
	//		if (flag) {
	//			continue;
	//		}
	//	}

	//	out.push_back(config);
	//	
	//}
	//
	return out;
}