#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>


const std::string testPath = "./data/";
const std::string testFiles[] = {
    "tsp_51_1",
    "tsp_70_1",
    "tsp_76_1",
    "tsp_76_2",
    "tsp_99_1",
    "tsp_100_1",
    "tsp_100_2",
    "tsp_100_3",
    "tsp_100_4",
    "tsp_100_5",
    "tsp_100_6",
};


struct Point {
    double x, y;
};


void calculate_total_way(const std::vector<Point>& points, std::vector<std::vector<double>>& distances) {
    for (size_t i = 0; i < points.size(); i++)
        for (size_t j = 0; j < points.size(); j++) {
            if (i != j)
                distances[i][j] = std::sqrt(
                    (points[i].x - points[j].x) * (points[i].x - points[j].x) +
                    (points[i].y - points[j].y) * (points[i].y - points[j].y));
        }
}

double calculate_way(const std::vector<std::vector<double>>& distances, const std::vector<size_t>& way) {
    double distance = 0.0;
    for (size_t i = 0; i < way.size(); i++)
        distance += distances[way[i]][way[(i + 1) % way.size()]];
    return distance;
}


double twp_opt(const std::vector<std::vector<double>>& distances, std::vector<size_t>& way) {
    double best_way;
    bool impr = true;
    while (impr) {
        impr = false;
        best_way = calculate_way(distances, way);
        for (size_t i = 1; i < way.size() - 1; ++i) {
            for (size_t j = i + 1; j < way.size(); ++j) {
                std::vector<size_t> new_way = way;
                std::reverse(new_way.begin() + i, new_way.begin() + j + 1);

                double new_dist = calculate_way(distances, new_way);
                if (new_dist < best_way) {
                    way = new_way;
                    best_way = new_dist;
                    impr = true;
                }
            }
        }
    }
    return best_way;
}


double three_opt(const std::vector<std::vector<double>>& distances, std::vector<size_t>& way) {
    size_t n = distances.size();
    double best_way;
    bool impr = true;
    while (impr) {
        impr = false;
        best_way = calculate_way(distances, way);
        for (size_t i = 0; i < n - 2; ++i) {
            for (size_t j = i + 1; j < n - 1; ++j) {
                for (size_t k = j + 1; k < n; ++k) {
                    std::vector<std::vector<size_t>> new_ways;
                    //перебираем возможные перестановки
                    new_ways.push_back({ way.begin(), way.begin() + i });
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - j), way.rbegin() + (n - i));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + j, way.end());

                    new_ways.push_back({ way.begin(), way.begin() + j });
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - k), way.rbegin() + (n - j));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + k, way.end());

                    new_ways.push_back({ way.rbegin() + (n - i), way.rend() });
                    new_ways.back().insert(new_ways.back().end(), way.rbegin(), way.rbegin() + (n - k));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + i, way.begin() + k);

                    new_ways.push_back({ way.begin(), way.begin() + i });
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - j), way.rbegin() + (n - i));
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - k), way.rbegin() + (n - j));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + k, way.end());

                    new_ways.push_back({ way.begin(), way.begin() + i });
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - k), way.rbegin() + (n - j));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + i, way.begin() + j);
                    new_ways.back().insert(new_ways.back().end(), way.begin() + k, way.end());

                    new_ways.push_back({ way.begin(), way.begin() + i });
                    new_ways.back().insert(new_ways.back().end(), way.begin() + j, way.begin() + k);
                    new_ways.back().insert(new_ways.back().end(), way.rbegin() + (n - j), way.rbegin() + (n - i));
                    new_ways.back().insert(new_ways.back().end(), way.begin() + k, way.end());

                    new_ways.push_back({ way.begin(), way.begin() + i });
                    new_ways.back().insert(new_ways.back().end(), way.begin() + j, way.begin() + k);
                    new_ways.back().insert(new_ways.back().end(), way.begin() + i, way.begin() + j);
                    new_ways.back().insert(new_ways.back().end(), way.begin() + k, way.end());

                    for (const auto& new_way : new_ways) {
                        double new_dist = calculate_way(distances, new_way);
                        if (new_dist < best_way) {
                            way = new_way; 
                            best_way = new_dist;
                            impr = true;
                        }
                    }
                }
            }
        }
    }
    return best_way;
}



int main() {


    std::cout << "TestFile" << "\t" << "\t" << "TSP Opt2 (distance/runtime)" << "\t" << "TSP Opt3 (distance/runtime)" << std::endl;

    for (std::string testFile : testFiles)
    {
        std::cout << testFile << "\t" << "\t";

        size_t n;
        std::vector<Point> points;

        std::ifstream in(testPath + testFile); // окрываем файл дл€ чтени€
        if (in.is_open()) {
            in >> n;
            double x, y;
            while (in >> x >> y)
                points.push_back(Point{ x, y });
            in.close();
        }
        else {
            std::cout << "File not found." << std::endl;
            return -1;
        }


        std::vector<std::vector<double>> distances(n, std::vector<double>(n, 0));
        calculate_total_way(points, distances);


        using std::chrono::high_resolution_clock;
        using std::chrono::duration_cast;
        using std::chrono::duration;
        using std::chrono::milliseconds;

        std::vector<size_t> way(n);

        auto t1 = high_resolution_clock::now();
        for (size_t i = 0; i < n; i++)    // инициализаци€ маршрута по умолчанию
            way[i] = i;
        double distance = twp_opt(distances, way);
        auto t2 = high_resolution_clock::now();
        auto ms = duration_cast<milliseconds>(t2 - t1);
        std::cout << distance << " / " << ms.count() << "ms" << "\t" << "\t" << "\t";


        t1 = high_resolution_clock::now();
        for (size_t i = 0; i < n; i++)    // инициализаци€ маршрута по умолчанию. ≈сли убрать, то будет использован маршрут после TSP 2opt
            way[i] = i;
        distance = three_opt(distances, way);
        t2 = high_resolution_clock::now();
        ms = duration_cast<milliseconds>(t2 - t1);
        std::cout << distance << " / " << ms.count() << "ms" << std::endl;
    }

    return 0;
}
