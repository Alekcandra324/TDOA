#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
//мне даны координаты точек D,F,E и разность ходов ( всего их 9 значений для каждой точки соот A,B,C)
// сначала я задам какие то координаты для неизвестных точек 
// мне понадобиться вычислить расстояние между точками, поэтому +1 функция, помимо основных.
// и основные функции реализуют вычисление ошибки и градиентный спуск.

// функция для вычисления расстояния между точками (например AD)- евклидова норма
double evclid(Vector2d point1, Vector2d point2) {
    return (point1 - point2).norm();
}

//вторая функция loss function L(A,B,C)
// для этого мне понадобяться предполагаемые координаты точек А,В и С , известные координаты D,F,E, и известные значения разности хода.
// использую метод наименьших квадратов
double loss_function(Vector2d A, Vector2d B, Vector2d C, const vector<Vector2d>& DEF_points, const vector<double>& distance_diff) {
    double mistake = 0.0;
    for (int i = 0; i < 3; i++) {
        Vector2d DEF = DEF_points[i];
        mistake += pow(evclid(A, DEF) - evclid(B, DEF) - distance_diff[i * 3], 2);
        mistake += pow(evclid(A, DEF) - evclid(C, DEF) - distance_diff[i * 3 + 1], 2);
        mistake += pow(evclid(B, DEF) - evclid(C, DEF) - distance_diff[i * 3 + 2], 2);
    }
    return mistake;

}

void Gradient_spusk(Vector2d& A, Vector2d& B, Vector2d& C, const vector<Vector2d>& DEF_points, const vector<double>& distance_diff, double step, int iterations) {
    double accuracy = 0.01;
    //предыдущий код работал неправильно, хотя математически мне казалось я реализовала его верно, однако ошибка с каждой итерацией росла.
    // необходимо было внимательнее читать википедию
    // буду реализовать Метод покоординатного спуска Гаусса — Зейделя
    //для этого я буду искать каждую координату поотдельности.
    for (int iter = 0; iter < iterations; iter++) {
        for (int coord = 0; coord < 2; coord++) { // 0 - x, 1 - y
            for (Vector2d* point : { &A, &B, &C }) {
                double orig_value = (*point)[coord];  // Запоминаем старое значение
                double grad = 0.0;
                // сложность, поняла из этой версии , что шаг итерации для каждой координаты должен быть отдельный и меняться.Иначе то в минус уходит , то начинает уходить от правильного ответа
                //буду думать как это объяснить машине
                for (int i = 0; i < 3; i++) {
                    Vector2d DEF = DEF_points[i];
                    double normA_DEF = (A - DEF).norm();
                    double normB_DEF = (B - DEF).norm();
                    double normC_DEF = (C - DEF).norm();

                    if (normA_DEF > 1e-8 && normB_DEF > 1e-8 && normC_DEF > 1e-8) {
                        if (point == &A) {
                            grad += 2 * ((normA_DEF - normB_DEF - distance_diff[i * 3]) * ((A - DEF)[coord] / normA_DEF) +
                                (normA_DEF - normC_DEF - distance_diff[i * 3 + 1]) * ((A - DEF)[coord] / normA_DEF));
                        }
                        if (point == &B) {
                            grad += 2 * ((normB_DEF - normC_DEF - distance_diff[i * 3 + 2]) * ((B - DEF)[coord] / normB_DEF) -
                                (normA_DEF - normB_DEF - distance_diff[i * 3]) * ((B - DEF)[coord] / normB_DEF));
                        }
                        if (point == &C) {
                            grad += 2 * (-(normA_DEF - normC_DEF - distance_diff[i * 3 + 1]) * ((C - DEF)[coord] / normC_DEF) -
                                (normB_DEF - normC_DEF - distance_diff[i * 3 + 2]) * ((C - DEF)[coord] / normC_DEF));
                        }
                    }
                }
                (*point)[coord] -= step * grad;
            }
        }
    }
}



// реализуем градиентный спуск
// для этого мне понадобиться вычислить значения частных производных L(A,B,C)  по каждой координате,
// функция должна выдать новые координаты , которые будут более правильными -> функция ошибки уменьшиться

int main() {
    

    // Известные точки D, E, F
    vector<Vector2d> DEF_points = { {40, 60}, {46, 52}, {34, 87} };

    // Разности хода сигнала 
    vector<double> distance_diff = { 5.78,20.89 , 15.11, 5.71, 19.75, 14.04, 5.70,21.45 , 15.75 };

    // Начальные догадки
    Vector2d A(10, 10), B(10, 20), C(20, 10);

    // Запускаем градиентный спуск
    Gradient_spusk(A, B, C, DEF_points, distance_diff, 0.0001, 10000);

    // Вывод результата
    cout << "A: (" << A.x() << ", " << A.y() << ")\n";
    cout << "B: (" << B.x() << ", " << B.y() << ")\n";
    cout << "C: (" << C.x() << ", " << C.y() << ")\n";
    return 0;
}