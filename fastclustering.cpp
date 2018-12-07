#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

namespace py = pybind11;


struct aminoAcid {
    std::string amino;
    std::vector<std::string> bases;
};

char getBaseName(std::set<char> baseSet) {
    if (baseSet == std::set<char>({'A'})) return 'A';
    if (baseSet == std::set<char>({'C'})) return 'C';
    if (baseSet == std::set<char>({'G'})) return 'G';
    if (baseSet == std::set<char>({'U'})) return 'U';
    if (baseSet == std::set<char>({'A', 'G'})) return 'R';
    if (baseSet == std::set<char>({'C', 'U'})) return 'Y';
    if (baseSet == std::set<char>({'G', 'C'})) return 'S';
    if (baseSet == std::set<char>({'A', 'U'})) return 'W';
    if (baseSet == std::set<char>({'G', 'U'})) return 'K';
    if (baseSet == std::set<char>({'A', 'C'})) return 'M';
    if (baseSet == std::set<char>({'C', 'G', 'U'})) return 'B';
    if (baseSet == std::set<char>({'A', 'G', 'U'})) return 'D';
    if (baseSet == std::set<char>({'A', 'C', 'U'})) return 'H';
    if (baseSet == std::set<char>({'A', 'C', 'G'})) return 'V';
    if (baseSet == std::set<char>({'A', 'C', 'G', 'U'})) return 'N';
    assert(false);
}

// https://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
template<typename T>
std::vector<T> sampleWithoutReplacement(std::vector<T> &array, int sampleSize) {
    std::vector<T> samples(sampleSize);

    int &n = sampleSize;
    int N = array.size();

    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    static std::random_device rnd;
    static std::mt19937 mt(rnd());
    static std::uniform_real_distribution<> uniform(0.0, 1.0); 

    while (m < n)
    {
        u = uniform(mt);

        if ( (N - t) * u >= n - m ) {
            t++;
        }
        else {
            samples[m] = array[t];
            t++; m++;
        }
    }
    return samples;
}

template<typename T>
T sampleOne(std::vector<T> &array) {
    static std::random_device rnd;
    static std::mt19937 mt(rnd());
    std::uniform_int_distribution<> uniform(0, array.size()-1); 
    return array[uniform(mt)];
}

int hammingDistance(const std::string &str1, const std::string &str2) {
    assert(str1.size() == str2.size());
    int dist = 0;
    for(size_t i=0; i<str1.size(); ++i) {
        if(str1[i] != str2[i]) dist++;
    }
    return dist;
}

std::tuple<int, std::string> minimumHammingDistance(const std::string &s, const std::vector<std::string> &a) {
    int minDist = std::numeric_limits<int>::max();
    std::string minStr;

    for(auto &t : a) {
        int d = hammingDistance(s, t);
        if (d < minDist) {
            minDist = d;
            minStr = t;
        }
    }
    return std::make_tuple(minDist, minStr);
}

// Assume array is sorted
template<typename T>
T mostCommonElement(const std::vector<T> &array) {
    auto max_val = *array.begin();
    auto test = max_val;
    int max_cnt = 0, cnt = 0;
    for(auto item : array) {
        if(item == test) {
            cnt++;
        }
        else {
            if(cnt > max_cnt) {
                max_val = test;
                max_cnt = cnt;
            }
            test = item;
            cnt = 1;
        }
    }
    return max_val;
}

std::vector<std::string> clustering(std::vector<aminoAcid> &aminos, int k) {
    std::random_device rnd;
    std::mt19937 mt(rnd()); 

    const int baseLen = aminos[0].bases[0].size();

    // Randomly nitialize centroids
    std::vector<aminoAcid> kAminos = sampleWithoutReplacement(aminos, k);
    std::vector<std::string> centroids;
    for(auto &amino : kAminos) {
        centroids.push_back(sampleOne(amino.bases));
    }

    std::vector<std::vector<std::string>> clusters(k);

    for (int step = 0; step < 5; ++step) {
        //std::cout << "Centroids:" << std::endl;
        //for (auto &c : centroids) {
        //    std::cout << c << std::endl;
        //}
        // assignment step
        clusters.clear(); clusters.resize(k);
        for (auto &amino : aminos) {
            int minDist = std::numeric_limits<int>::max();
            int minCluster = -1;
            std::string minBase;
            for (size_t c = 0; c < centroids.size(); ++c) {
                auto t = minimumHammingDistance(centroids[c], amino.bases);
                int d = std::get<0>(t);
                std::string b = std::get<1>(t);
                if (d < minDist) {
                    minBase = b;
                    minCluster = c;
                    minDist = d;
                }
            }
            clusters[minCluster].push_back(minBase);
        }

        //for(auto &cluster : clusters) {
        //    std::cout << "cluster size: " << cluster.size() << std::endl;
        //    for (auto &c : cluster) {
        //        std::cout << c << ",";
        //    }
        //    std::cout << std::endl;
        //}
        
        // Update step
        centroids.clear(); centroids.resize(k);
        for (size_t i=0; i<clusters.size(); ++i) {
            auto &cluster = clusters[i];
            assert(!cluster.empty());
            for (int j=0; j<baseLen; ++j) {
                std::vector<char> tmp(cluster.size());
                for (size_t p=0; p<cluster.size(); ++p) {
                    tmp[p] = cluster[p][j];
                }
                std::sort(tmp.begin(), tmp.end());
                centroids[i] += mostCommonElement(tmp);
            }
        }
    }

    std::vector<std::string> primers(k);
    for (size_t i=0; i<clusters.size(); ++i) {
        auto &cluster = clusters[i];
        assert(!cluster.empty());
        for (int j=0; j<baseLen; ++j) {
            std::set<char> B;
            for (size_t p=0; p<cluster.size(); ++p) {
                B.insert(cluster[p][j]);
            }
            primers[i] += getBaseName(B);
        }
    }

    return primers;
}

PYBIND11_MODULE(fastclustering, m) {
    m.doc() = "Fast clustering module";
    py::class_<aminoAcid>(m, "aminoAcid")
        .def(py::init<>())
        .def_readwrite("amino", &aminoAcid::amino)
        .def_readwrite("bases", &aminoAcid::bases);
    m.def("clustering", &clustering);
}
