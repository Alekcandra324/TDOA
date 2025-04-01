/*#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <unordered_map>
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

void Gradient_spusk(Vector2d& A, Vector2d& B, Vector2d& C, const vector<Vector2d>& DEF_points, const vector<double>& distance_diff,  int iterations) {
    //double accuracy = 0.001;
    unordered_map<Vector2d*, array<double, 2>> step_sizes = {
    { &A, {0.01, 0.01} },
    { &B, {0.001, 0.001} },
    { &C, {0.001, 0.001} }//добавила отдельные шаги

    };

    //предыдущий код работал неправильно, хотя математически мне казалось я реализовала его верно, однако ошибка с каждой итерацией росла.
    // необходимо было внимательнее читать википедию
    // буду реализовать Метод покоординатного спуска Гаусса — Зейделя
    //для этого я буду искать каждую координату поотдельности.
    for (int iter = 0; iter < iterations; iter++) {
        double max_change = 0.0;
        

        for (int coord = 0; coord < 2; coord++) { // 0 - x, 1 - y
            for (Vector2d* point : { &C, &B,&A }) {

                double grad = 0.0;


                // сложность, поняла из этой версии , что шаг итерации для каждой координаты должен быть отдельный и меняться.Иначе то в минус уходит , то начинает уходить от правильного ответа

                for (int i = 0; i < 3; i++) {
                    Vector2d DEF = DEF_points[i];
                    double normA_DEF = (A - DEF).norm();
                    double normB_DEF = (B - DEF).norm();
                    double normC_DEF = (C - DEF).norm();
                    if (point == &A) {
                        grad += 2 * ((normA_DEF - normB_DEF - distance_diff[i * 3]) * ((A - DEF)[coord] / normA_DEF) +
                            (normA_DEF - normC_DEF - distance_diff[i * 3 + 1]) * ((A - DEF)[coord] / normA_DEF));
                    }
                    if (point == &B) {
                        grad += 2 * ((normB_DEF - normC_DEF - distance_diff[i * 3 + 2]) * ((B - DEF)[coord] / normB_DEF) -
                            (normA_DEF - normB_DEF - distance_diff[i * 3]) * ((B - DEF)[coord] / normB_DEF));
                    }
                    if (point == &C) {
                        grad += 2 * ((normA_DEF - normC_DEF - distance_diff[i * 3 + 1]) * ((C - DEF)[coord] / normC_DEF) -
                            (normB_DEF - normC_DEF - distance_diff[i * 3 + 2]) * ((C - DEF)[coord] / normC_DEF));
                    }

                }
                double orig = (*point)[coord];
                double cost = loss_function(A, B, C, DEF_points, distance_diff);
                (*point)[coord] -= step_sizes[point][coord] * grad;
                double new_cost = loss_function(A, B, C, DEF_points, distance_diff);
                if (new_cost > cost) {
                    step_sizes[point][coord] *= 0.1;
                    (*point)[coord] = orig;
                    (*point)[coord] -= step_sizes[point][coord] * grad;
                }
                
                cout << "Updated " << ((point == &A) ? "A" : (point == &B) ? "B" : "C")
                    << "[" << coord << "] to " << (*point)[coord] << ", new loss: " << new_cost << ", past loss: " << cost <<"grad "<<grad<< " step size " << step_sizes[point][coord] << "\n";         

            }
        }
        cout << "Iter " << iter << " Loss: " << loss_function(A, B, C, DEF_points, distance_diff) << "\n";

       // if (max_change < accuracy) break; // Условие остановки
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
    Vector2d A(3, 3), B(8, 8), C(10, 19);

    // Запускаем градиентный спуск
    Gradient_spusk(A, B, C, DEF_points, distance_diff,  500);

    // Вывод результата
    cout << "A: (" << A.x() << ", " << A.y() << ")\n";
    cout << "B: (" << B.x() << ", " << B.y() << ")\n";
    cout << "C: (" << C.x() << ", " << C.y() << ")\n";
    return 0;
    
}
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <unordered_map>
using namespace std;
using namespace Eigen;
//изучив немного теории, решила начать сначала, написать программу для нахождения координат одной антенны, затем реализовать для трех.
//интеренет говорит, что для нелинейной системы лучше подходит или обычный GD  или SGD, а я похоже запуталась, поэтому реализую просто GD


//первые две функции не отличаются от предыдущего кода, почти.Разница в том, что гиперболу я рассчитываю как разницу между расстояниями (от А до D) и (от А до Е) например, а не как разницу расстояний двух антенн до одного и того же источника.
double evclid(Vector2d point1, Vector2d point2) {
    return (point1 - point2).norm();
}

double loss_function(Vector2d A, Vector2d D, Vector2d E, Vector2d F, const vector<double>& distance_diff) {
    double mistake = 0.0;
    //0- D,1-E,2-F
    mistake += pow(evclid(A, D) - evclid(A, E) - distance_diff[0], 2);
    mistake += pow(evclid(A, D) - evclid(A, F) - distance_diff[1], 2);
    mistake += pow(evclid(A, E) - evclid(A, F) - distance_diff[2], 2);
    return mistake;
}



void Gradient_spusk(Vector2d& A, Vector2d D, Vector2d E, Vector2d F, const vector<double>& distance_diff, int iterations) {
    double step_sizes[2] = { 0.001, 0.0001 };

    //θ(t + 1) = θ(t)−α⋅∇f(θ(t)) реализовать надо вот это
    for (int iter = 0; iter < iterations; iter++) {
        //получается  teta  это или x  или y 
        for (int coord = 0; coord < 2; coord++) {
            // тогда выходит три градиента за один проход для  x и y. 
            double old_loss = loss_function(A, D, E, F, distance_diff);

            double grad = 0.0;
            double normA_D = (A - D).norm();
            double normA_E = (A - E).norm();
            double normA_F = (A - F).norm();
            grad += 2 * (normA_D - normA_E - distance_diff[0]) * ((A - D)[coord] / normA_D - (A - E)[coord] / normA_E);
            grad += 2 * (normA_D - normA_F - distance_diff[1]) * ((A - D)[coord] / normA_D - (A - F)[coord] / normA_F);
            grad += 2 * (normA_E - normA_F - distance_diff[2]) * ((A - E)[coord] / normA_E - (A - F)[coord] / normA_F);
            double old_coord = A[coord];
            A[coord] -= step_sizes[coord] * grad;
           
        

        }
    }

}

int main() {

    // Разности хода сигнала 
    vector<double> distance_diff;
    // Известные точки D, E, F
    Vector2d D(5, 1), E(5, 8), F(7, 9);
    //Известные координаты для А
    Vector2d A_v(10, 6);
    distance_diff.push_back((A_v - D).norm() - (A_v - E).norm());
    distance_diff.push_back((A_v - D).norm() - (A_v - F).norm());
    distance_diff.push_back((A_v - E).norm() - (A_v - F).norm());

    // Начальные догадки
    Vector2d A(3, 7);
    // нашел просто идеально координаты A_v(5, 1) при шаге для каждой 0.00001 и 10000 итераций и ошибка = 10^(-10)
    // при тех же параметрах для A_v(4, 6) плохо работает.
    // Запускаем градиентный спуск
    Gradient_spusk(A, D, E, F, distance_diff, 100000);

    // Вывод результата
    cout << "A: (" << A.x() << ", " << A.y() << ")\n";
    cout << loss_function(A, D, E, F, distance_diff);
    return 0;

}