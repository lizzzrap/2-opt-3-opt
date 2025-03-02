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


void calculateTotalDistances(const std::vector<Point>& points, std::vector<std::vector<double>>& distances) {
    for (size_t i = 0; i < points.size(); i++) {
        for (size_t j = 0; j < points.size(); j++) {
            if (i != j)
                distances[i][j] = std::sqrt(
                    (points[i].x - points[j].x) * (points[i].x - points[j].x) +
                    (points[i].y - points[j].y) * (points[i].y - points[j].y));
        }
    }
}

double calculateRouteDistance(const std::vector<std::vector<double>>& distances, const std::vector<size_t>& route) {
    double distance = 0.0;
    for (size_t i = 0; i < route.size(); i++) {
        distance += distances[route[i]][route[(i + 1) % route.size()]];
    }
    return distance;
}


double twoOpt(const std::vector<std::vector<double>>& distances, std::vector<size_t>& route) {
    double bestDistance;
    bool improvement = true;
    while (improvement) {
        improvement = false;
        bestDistance = calculateRouteDistance(distances, route);
        for (size_t i = 1; i < route.size() - 1; ++i) {
            for (size_t j = i + 1; j < route.size(); ++j) {
                // Создаем новую последовательность, меняя порядок между i и j
                std::vector<size_t> newRoute = route;
                std::reverse(newRoute.begin() + i, newRoute.begin() + j + 1);

                double newDistance = calculateRouteDistance(distances, newRoute);
                if (newDistance < bestDistance) {
                    route = newRoute; // Обновляем точки
                    bestDistance = newDistance;
                    improvement = true; // Найдено улучшение
                }
            }
        }
    }
    return bestDistance;
}


double threeOpt(const std::vector<std::vector<double>>& distances, std::vector<size_t>& route) {
    size_t n = distances.size();
    double bestDistance;
    bool improvement = true;
    while (improvement) {
        improvement = false;
        bestDistance = calculateRouteDistance(distances, route);
        for (size_t i = 0; i < n - 2; ++i) {
            for (size_t j = i + 1; j < n - 1; ++j) {
                for (size_t k = j + 1; k < n; ++k) {
                    // Создаем новые маршруты
                    std::vector<std::vector<size_t>> newRoutes;

                    // 1. Перевернуть порядок между i и j, оставить k
                    newRoutes.push_back({ route.begin(), route.begin() + i });
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - j), route.rbegin() + (n - i));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + j, route.end());

                    // 2. Перевернуть порядок между j и k, оставить i
                    newRoutes.push_back({ route.begin(), route.begin() + j });
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - k), route.rbegin() + (n - j));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + k, route.end());

                    // 3. Перевернуть порядок между i и k, оставить j
                    newRoutes.push_back({ route.rbegin() + (n - i), route.rend() });
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin(), route.rbegin() + (n - k));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + i, route.begin() + k);

                    // 4. k->i + j->i + k->j
                    newRoutes.push_back({ route.begin(), route.begin() + i });
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - j), route.rbegin() + (n - i));
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - k), route.rbegin() + (n - j));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + k, route.end());

                    // 5.  k->i + k->j + i->j
                    newRoutes.push_back({ route.begin(), route.begin() + i });
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - k), route.rbegin() + (n - j));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + i, route.begin() + j);
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + k, route.end());

                    // 6.  k->i + j->k + j->i
                    newRoutes.push_back({ route.begin(), route.begin() + i });
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + j, route.begin() + k);
                    newRoutes.back().insert(newRoutes.back().end(), route.rbegin() + (n - j), route.rbegin() + (n - i));
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + k, route.end());

                    // 7.  k->i + j->k + i->j
                    newRoutes.push_back({ route.begin(), route.begin() + i });
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + j, route.begin() + k);
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + i, route.begin() + j);
                    newRoutes.back().insert(newRoutes.back().end(), route.begin() + k, route.end());

                    // Проверка новых маршрутов
                    for (const auto& newRoute : newRoutes) {
                        double newDistance = calculateRouteDistance(distances, newRoute);
                        if (newDistance < bestDistance) {
                            route = newRoute; // Обновляем точки
                            bestDistance = newDistance;
                            improvement = true; // Найдено улучшение
                        }
                    }
                }
            }
        }
    }
    return bestDistance;
}



int main() {


    std::cout << "TestFile" << "\t" << "\t" << "TSP Opt2 (distance/runtime)" << "\t" << "TSP Opt3 (distance/runtime)" << std::endl;

    for (std::string testFile : testFiles)
    {
        std::cout << testFile << "\t" << "\t";

        size_t n;
        std::vector<Point> points;

        std::ifstream in(testPath + testFile); // окрываем файл для чтения
        if (in.is_open()) {
            in >> n;
            double x, y;
            while (in >> x >> y) {
                points.push_back(Point{ x, y });
            }
            in.close();
        }
        else {
            std::cout << "File not found." << std::endl;
            return -1;
        }


        std::vector<std::vector<double>> distances(n, std::vector<double>(n, 0));
        calculateTotalDistances(points, distances);


        using std::chrono::high_resolution_clock;
        using std::chrono::duration_cast;
        using std::chrono::duration;
        using std::chrono::milliseconds;

        std::vector<size_t> route(n);

        auto t1 = high_resolution_clock::now();
        for (size_t i = 0; i < n; i++) {    // инициализация маршрута по умолчанию
            route[i] = i;
        }
        double distance = twoOpt(distances, route);
        auto t2 = high_resolution_clock::now();
        auto ms = duration_cast<milliseconds>(t2 - t1);
        std::cout << distance << " / " << ms.count() << "ms" << "\t" << "\t" << "\t";


        t1 = high_resolution_clock::now();
        for (size_t i = 0; i < n; i++) {    // инициализация маршрута по умолчанию. Если убрать, то будет использован маршрут после TSP 2opt
            route[i] = i;
        }
        distance = threeOpt(distances, route);
        t2 = high_resolution_clock::now();
        ms = duration_cast<milliseconds>(t2 - t1);
        std::cout << distance << " / " << ms.count() << "ms" << std::endl;
    }

    return 0;
}
