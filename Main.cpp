#include <bits/stdc++.h>
using namespace std;

void clearScreen()
{
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

void JacobiIteration()
{
    // 27 6 -1 85
    // 6 15 2 72
    // 1 2 54 110

    cout << "enter size" << endl;
    int n;
    cin >> n;
    vector<vector<double>> AugmentedMatrix(n, vector<double>(n + 1));
    vector<double> x(n, 0.0);

    cout << "Enter the augmented matrix (coefficients and constants):" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            cin >> AugmentedMatrix[i][j];
        }
    }
    double tolerance = 1e-6;
    int maxIterations = 1000;

    clearScreen();

    vector<double> x_new(n, 0.0);

    for (int iter = 1; iter <= maxIterations; ++iter)
    {
        bool converged = true;

        for (int i = 0; i < n; ++i)
        {
            x_new[i] = AugmentedMatrix[i][n];
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    x_new[i] -= AugmentedMatrix[i][j] * x[j];
                }
            }
            x_new[i] /= AugmentedMatrix[i][i];
            if (abs(x_new[i] - x[i]) > tolerance)
            {
                converged = false;
            }
        }

        x = x_new;

        if (converged)
        {
            cout << "Converged after " << iter << " iterations." << endl;
            break;
        }
    }

    cout << "Final solution: ";
    for (int i = 0; i < n; ++i)
    {
        cout << fixed << setprecision(6) << x[i] << " ";
    }
    cout << endl;

    // return 0;
}

// void GaussSeidelIteration()
// {
//     // Example augmented matrix:
//     // 3 2 -1 11
//     // 2 -3 1 7
//     // 5 1 -2 12
//     cout << "Enter size" << endl;
//     int n;
//     cin >> n;
//     vector<vector<double>> AugmentedMatrix(n, vector<double>(n + 1));
//     vector<double> x(n, 0.0);

//     cout << "Enter the augmented matrix (coefficients and constants):" << endl;
//     for (int i = 0; i < n; ++i)
//     {
//         for (int j = 0; j <= n; ++j)
//         {
//             cin >> AugmentedMatrix[i][j];
//         }
//     }
//     clearScreen();

//     double tolerance = 1e-6;
//     int maxIterations = 100000;

//     for (int iter = 1; iter <= maxIterations; ++iter)
//     {
//         bool converged = true;

//         for (int i = 0; i < n; ++i)
//         {
//             double sum = AugmentedMatrix[i][n];
//             for (int j = 0; j < n; ++j)
//             {
//                 if (j < i)
//                 {
//                     sum -= AugmentedMatrix[i][j] * x[j];
//                 }
//                 else
//                 {
//                     sum -= AugmentedMatrix[i][j] * x[j];
//                 }
//             }
//             double x_new = sum / AugmentedMatrix[i][i];

//             if (abs(x_new - x[i]) > tolerance)
//             {
//                 converged = false;
//             }

//             x[i] = x_new;
//         }

//         if (converged)
//         {
//             cout << "Converged after " << iter << " iterations." << endl;
//             break;
//         }
//     }

//     cout << "Final solution: ";
//     for (int i = 0; i < n; ++i)
//     {
//         cout << fixed << setprecision(6) << x[i] << " ";
//     }
//     cout << endl;
// }

void GaussSeidelIteration()
{
    // Example augmented matrix:
    // 3 2 -1 11
    // 2 -3 1 7
    // 5 1 -2 12
    cout << "Enter size" << endl;
    int n;
    cin >> n;
    vector<vector<double>> AugmentedMatrix(n, vector<double>(n + 1));
    vector<double> x(n, 0.0);

    cout << "Enter the augmented matrix (coefficients and constants):" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            cin >> AugmentedMatrix[i][j];
        }
    }
    clearScreen();

    double tolerance = 1e-6;
    int maxIterations = 100000;

    for (int iter = 1; iter <= maxIterations; ++iter)
    {
        bool converged = true;

        for (int i = 0; i < n; ++i)
        {
            double sum = AugmentedMatrix[i][n];
            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    sum -= AugmentedMatrix[i][j] * x[j];
                }
            }
            double x_new = sum / AugmentedMatrix[i][i];

            if (abs(x_new - x[i]) > tolerance)
            {
                converged = false;
            }

            x[i] = x_new;
        }

        if (converged)
        {
            cout << "Converged after " << iter << " iterations." << endl;
            break;
        }
    }

    cout << "Final solution: ";
    for (int i = 0; i < n; ++i)
    {
        cout << fixed << setprecision(6) << x[i] << " ";
    }
    cout << endl;
}

void GaussElimination()
{
    cout << "Enter size" << endl;
    int n;
    cin >> n;

    vector<vector<float>> arr(n, vector<float>(n + 1));
    cout << "Enter augmented matrix" << endl;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cin >> arr[i][j];
        }
    }

    clearScreen();

    vector<float> x(n);
    float ratio;

    for (int i = 0; i < n - 1; i++)
    {
        if (arr[i][i] == 0)
        {
            cerr << "Error: Division by zero, zero pivot encountered." << endl;
            return;
        }

        for (int j = i + 1; j < n; j++)
        {
            ratio = arr[j][i] / arr[i][i];
            for (int k = 0; k <= n; k++)
            {
                arr[j][k] -= arr[i][k] * ratio;
            }
        }
    }

    x[n - 1] = arr[n - 1][n] / arr[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = arr[i][n];
        for (int j = i + 1; j < n; j++)
        {
            x[i] -= arr[i][j] * x[j];
        }
        x[i] /= arr[i][i];
    }

    cout << "Solution: ";

    for (auto it : x)
    {
        cout << it << " ";
    }
    cout << endl;
}

void GaussJordanElimination()
{
    cout << "Enter size" << endl;
    int n;
    cin >> n;

    vector<vector<float>> arr(n, vector<float>(n + 1));
    cout << "Enter augmented matrix" << endl;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cin >> arr[i][j];
        }
    }
    clearScreen();

    vector<float> x(n);
    float ratio;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            ratio = arr[j][i] / arr[i][i];
            for (int k = 0; k <= n; k++)
            {
                arr[j][k] -= arr[i][k] * ratio;
            }
        }
    }

    for (int i = n - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            ratio = arr[j][i] / arr[i][i];
            for (int k = 0; k <= n; k++)
            {
                arr[j][k] -= arr[i][k] * ratio;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        x[i] = arr[i][n] / arr[i][i];
    }

    for (auto it : x)
    {
        cout << it << " ";
    }
}

double fx(double x, vector<int> &arr)
{
    double ans = 0;
    int n = arr.size();
    for (int i = 0; i < n; i++)
    {
        ans += arr[i] * pow(x, n - 1 - i);
    }
    return ans;
}

void UXsubstitution(vector<vector<float>> &U, vector<float> &X, vector<float> &Y)
{
    int n = U.size();

    X[n - 1] = Y[n - 1] / U[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++)
        {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}

void LYsubstitution(vector<vector<float>> &L, vector<float> &Y, vector<float> &B)
{
    int n = L.size();

    Y[0] = B[0] / L[0][0];

    for (int i = 1; i <= n - 1; i++)
    {
        Y[i] = B[i];
        for (int j = i - 1; j >= 0; j--)
        {
            Y[i] -= L[i][j] * Y[j];
        }
        Y[i] /= L[i][i];
    }
}

void decompose(vector<vector<float>> &L, vector<vector<float>> &U)
{
    int n = U.size();

    // making L and U
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            float diag = U[j][i] / U[i][i];
            L[j][i] = diag;
            for (int k = i; k < n; k++)
            {
                U[j][k] -= U[i][k] * diag;
            }
        }
    }
}

void LU()
{
    int n;
    cout << "enter size" << endl;
    cin >> n;
    vector<vector<float>> arr(n, vector<float>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> arr[i][j];
        }
    }
    vector<float> B(n);
    for (int i = 0; i < n; i++)
    {
        cin >> B[i];
    }

    clearScreen();

    vector<vector<float>> U = arr;
    vector<vector<float>> L(n, vector<float>(n, 0));

    for (int i = 0; i < n; i++)
    {
        L[i][i] = 1;
    }

    vector<float> Y(n);
    vector<float> X(n);

    decompose(L, U);
    LYsubstitution(L, Y, B);
    UXsubstitution(U, X, Y);

    for (auto it : X)
    {
        cout << it << " ";
    }
}

void Inverse(vector<vector<float>> &arr, vector<vector<float>> &unit, int n)
{
    vector<float> x(n);
    float ratio;

    for (int i = 0; i < n - 1; i++)
    {
        if (arr[i][i] == 0)
        {
            cerr << "Error: Division by zero, zero pivot encountered." << endl;
            return;
        }

        for (int j = i + 1; j < n; j++)
        {
            ratio = arr[j][i] / arr[i][i];
            for (int k = 0; k < n; k++)
            {
                arr[j][k] -= arr[i][k] * ratio;
                unit[j][k] -= unit[i][k] * ratio;
            }
        }
    }
    for (int i = n - 1; i > 0; i--)
    {
        if (arr[i][i] == 0)
        {
            cerr << "Error: Division by zero, zero pivot encountered." << endl;
            return;
        }

        for (int j = i - 1; j >= 0; j--)
        {
            ratio = arr[j][i] / arr[i][i];
            for (int k = 0; k < n; k++)
            {
                arr[j][k] -= arr[i][k] * ratio;
                unit[j][k] -= unit[i][k] * ratio;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        float diag = arr[i][i];
        for (int j = 0; j < n; j++)
        {
            arr[i][j] /= diag;
            unit[i][j] /= diag;
        }
    }

    for (auto it : unit)
    {
        for (int i = 0; i < it.size(); i++)
        {
            cout << it[i] << " ";
        }
        cout << endl;
    }
}
void InverseMatrix()
{
    int n;
    cout << "Enter Size" << endl;
    cin >> n;
    cout << "enter matrix" << endl;
    vector<vector<float>> arr(n, vector<float>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> arr[i][j];
        }
    }

    vector<vector<float>> unit(n, vector<float>(n, 0));
    for (int i = 0; i < n; i++)
    {
        unit[i][i] = 1;
    }

    Inverse(arr, unit, n);
}

void Bisection()
{
    int p;
    cout << "Enter power: ";
    cin >> p;

    vector<int> arr(p + 1);
    cout << "Enter coefficients: ";
    for (int i = 0; i < arr.size(); i++)
    {
        cin >> arr[i];
    }

    double xMax = ceil(sqrt(pow(arr[1] / arr[0], 2) - 2 * (arr[2] / arr[0])));
    double r1 = -xMax;
    double r2 = xMax;

    double a = r1, b = r2;
    bool intervalFound = false;

    for (int i = r1; i <= r2; i++)
    {
        if (fx(i, arr) == 0)
        {
            cout << "root = " << i << endl;
            return;
        }
        else if ((fx(i, arr) * fx(i + 1, arr)) < 0)
        {
            a = i;
            b = i + 1;
            intervalFound = true;
            break;
        }
    }

    double c;
    double tolerance = 1e-6;
    int max_iterations = 2000;
    int iteration = 0;

    while ((b - a) >= tolerance && iteration < max_iterations)
    {
        c = (a + b) / 2;

        if (fabs(fx(c, arr)) < tolerance)
        {
            cout << "Exact root is found within tolerance." << endl;
            cout << "root = " << c << endl;
            return;
        }

        if (fx(a, arr) * fx(c, arr) < 0)
            b = c;
        else
            a = c;

        iteration++;
    }

    cout << "root = " << c << endl;
}

void FalsePosition()
{
    int p;
    int t_roots = 0;
    vector<float> roots;
    cout << "Enter power: ";
    cin >> p;

    vector<int> arr(p + 1);
    cout << "Enter coefficients: ";
    for (int i = 0; i < arr.size(); i++)
    {
        cin >> arr[i];
    }

    double xMax = ceil(sqrt(pow(arr[1] / (double)arr[0], 2) - 2 * (arr[2] / (double)arr[0])));
    double r1 = -xMax;
    double r2 = xMax;

    double a = r1, b = r2;
    bool intervalFound = false;

    for (int i = r1; i <= r2; i++)
    {
        if (fx(i, arr) == 0)
        {
            cout << "root is: " << i << endl;
            break;
        }

        else if ((fx(i, arr) * fx(i + 1, arr)) < 0)
        {
            a = i;
            b = i + 1;
            intervalFound = true;
            break;
        }
    }
    if (intervalFound)
    {

        double tolerance = 1e-6;
        int max_iterations = 1000;
        double c, fa = fx(a, arr), fb = fx(b, arr);
        int prev_c = max_iterations;
        // for (int i = 1; abs(b - a) >= tolerance && i <= max_iterations; i++)
        while (1)
        {
            c = (a * fb - b * fa) / (fb - fa);
            double fc = fx(c, arr);

            if (fabs(fx(c, arr)) < tolerance)
            {
                cout << "Exact root is found within tolerance." << endl;
                cout << "root = " << c << endl;
                return;
            }

            else if (abs(prev_c - c) < tolerance)
            {
                cout << "Exact root is found within tolerance." << endl;
                cout << "root = " << c << endl;
                return;
            }
            else if (fa * fc < 0)
            {
                b = c;
            }
            else if (fa * fc > 0)
                a = c;
            prev_c = c;
        }

        cout << "root = " << c << endl;
    }
}

double fx1(double a, vector<double> &arr)
{
    int n = arr.size();
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += arr[i] * pow(a, n - 1 - i);
    }
    return sum;
}
double diff(double a, vector<double> &arr)
{
    int n = arr.size();
    double sum = 0;
    for (int i = 0; i < n - 1; i++)
    {
        sum += arr[i] * (n - 1 - i) * pow(a, n - i - 2);
    }
    return sum;
}

double newtonRaphson(vector<double> &arr, double a)
{
    double tolerance = 1e-6;
    int maxIter = 1000;
    double nextX = 0;
    double x = a;

    for (int i = 0; i < maxIter; i++)
    {

        nextX = x - fx1(x, arr) / diff(x, arr);

        if (abs(x - nextX) < tolerance)
        {
            return nextX;
        }

        x = nextX;
    }
    return x;
}

void deflate(vector<double> &arr, double root)
{
    int n = arr.size();
    vector<double> newArr(n - 1);
    newArr[0] = arr[0];
    for (int i = 1; i < n - 1; i++)
    {
        newArr[i] = arr[i] + newArr[i - 1] * root;
    }
    arr = newArr;
}

void Newtonraphson()
{
    int p;
    cout << "Enter Power" << endl;
    cin >> p;
    vector<double> arr(p + 1);

    for (int i = 0; i <= p; i++)
    {
        cin >> arr[i];
    }

    clearScreen();

    vector<double> roots;
    double x = 10;
    while (p > 0)
    {
        double root = newtonRaphson(arr, x);
        roots.push_back(root);
        deflate(arr, root);
        p--;
    }

    for (auto it : roots)
    {
        cout << it << " ";
    }
}

double fx(double a, vector<double> &arr)
{
    double sum = 0;
    int n = arr.size();

    for (int i = 0; i < n; i++)
    {
        sum += arr[i] * pow(a, n - 1 - i);
    }
    return sum;
}

double secant(vector<double> &arr, double a, double b)
{
    double tolerance = 1e-6;
    int maxIter = 1000;
    double x1 = a, x2 = b;
    double x3 = 0;

    for (int i = 0; i < maxIter; i++)
    {

        x3 = (x1 * fx(x2, arr) - x2 * fx(x1, arr)) / (fx(x2, arr) - fx(x1, arr));
        if (abs(x1 - x2) < tolerance)
        {
            return x3;
        }

        x1 = x2;
        x2 = x3;
    }
    return x2;
}

void Secant()
{
    int p;
    cout << "Enter Power" << endl;
    cin >> p;
    vector<double> arr(p + 1);
    for (int i = 0; i < p + 1; i++)
    {
        cin >> arr[i];
    }

    clearScreen();

    vector<double> roots;
    double x1 = 10, x2 = 12;
    while (p--)
    {
        double root = secant(arr, x1, x2);
        roots.push_back(root);
        deflate(arr, root);
    }
    for (auto it : roots)
    {
        cout << it << " ";
    }
}

double f(double x, double y)
{
    return x + y; // Example: dy/dx = x + y
}
void RungeKutta()
{
    double x0, y0, h, x_target;

    cout << "Enter initial values for x and y: ";
    cin >> x0 >> y0;
    cout << "Enter step size h: ";
    cin >> h;
    cout << "Enter target x value: ";
    cin >> x_target;

    clearScreen();

    double x = x0;
    double y = y0;

    cout << "x" << setw(10) << "y" << endl;
    cout << fixed << setprecision(6);

    while (x < x_target)
    {
        double k1 = h * f(x, y);
        double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
        double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
        double k4 = h * f(x + h, y + k3);

        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        x = x + h;

        cout << x << setw(10) << y << endl;
    }
}

int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

double newtonForward(vector<double> &x, vector<double> &y, double value)
{
    int n = x.size();
    vector<vector<double>> diff(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
        }
    }

    double result = y[0];
    double u = (value - x[0]) / (x[1] - x[0]);
    double term = u;

    for (int i = 1; i < n; i++)
    {
        result += (term * diff[0][i]) / fact(i);
        term *= (u - i);
    }

    return result;
}

double newtonBackward(vector<double> &x, vector<double> &y, double value)
{
    int n = x.size();
    vector<vector<double>> diff(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++)
    {
        for (int i = n - 1; i >= j; i--)
        {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }

    double result = y[n - 1];
    double u = (value - x[n - 1]) / (x[1] - x[0]);
    double term = 1.0;

    for (int i = 1; i < n; i++)
    {
        term *= (u + (i - 1));
        result += (term * diff[n - 1][i]) / fact(i);
    }

    return result;
}

void NewtonForwardBackwardInterpolation()
{
    int n;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter data points (x and y values):\n";
    for (int i = 0; i < n; i++)
    {
        cin >> x[i] >> y[i];
    }

    double value;
    cout << "Enter the value of x for which you want to interpolate: ";
    cin >> value;

    clearScreen();

    double resultForward = newtonForward(x, y, value);
    cout << "The interpolated value at x = " << value << " using Forward Interpolation is: " << fixed << setprecision(5) << resultForward << endl;

    double resultBackward = newtonBackward(x, y, value);
    cout << "The interpolated value at x = " << value << " using Backward Interpolation is: " << fixed << setprecision(5) << resultBackward << endl;
}

int main()
{
    int p;
    while (1)
    {
        cout << "Press 0 to Inverse a matrix" << endl;
        cout << "Press 1 to solve system of linear equations" << endl;
        cout << "Press 2 to solve polynomial equation" << endl;
        cout << "Press 3 to use Runge Kutta Method" << endl;
        cout << "Press 4 to use Newton Forward Backward Interpolation method" << endl;
        cout << "press 5 to terminate the process" << endl;
        cin >> p;

        clearScreen();
        if (p == 0)
        {
            // 3
            //  2 0 0
            //  0 3 0
            //  0 0 3
            InverseMatrix();
            return 0;
        }
        if (p == 1)
        {
            int p1;

            cout << "Press 1 to use Jacobi Iterative Method" << endl;
            cout << "Press 2 to use Gauss Seidel Iterative Method" << endl;
            cout << "Press 3 to use Gauss Elimination Method" << endl;
            cout << "Press 4 to use Gauss Jordan Elimination method" << endl;
            cout << "Press 5 to use LU Factorization" << endl;
            cout << "Press 6 to terminate the process" << endl;
            cin >> p1;

            clearScreen();

            if (p1 == 1)
            {
                JacobiIteration();
                cout << endl;
                return 0;
                // 3 2 -1 11
                // 2 -3 1 7
                // 5 1 -2 12
            }
            else if (p1 == 2)
            {
                GaussSeidelIteration();
                cout << endl;
                return 0;
                // 3 2 -1 11
                // 2 -3 1 7
                // 5 1 -2 12
            }
            else if (p1 == 3)
            {
                GaussElimination();
                cout << endl;
                return 0;

                // 3
                // 1 -1 2 3
                // 1 2 3 5
                // 3 -4 -5 -13
            }
            else if (p1 == 4)
            {
                GaussJordanElimination();
                cout << endl;
                return 0;

                // 3
                // 1 -1 2 3
                // 1 2 3 5
                // 3 -4 -5 -13
            }
            else if (p1 == 5)
            {
                LU();
                cout << endl;
                return 0;

                // 3
                // 1 5 1
                // 2 1 3
                // 3 1 4

                // 14 13 17
            }
            else if (p1 == 6)
            {
                cout << "The Process Is Terminated" << endl;
                return 0;
            }
            else
            {
                cout << "Invalid Input" << endl;
                return 0;
            }
        }

        else if (p == 2)
        {
            int p2;

            cout << "Press 1 to use Bisection Method" << endl;
            cout << "Press 2 to use False Position Method" << endl;
            cout << "Press 3 to use Newton Raphson Method" << endl;
            cout << "Press 4 to use Secant Method" << endl;
            cout << "Press 5 to terminate the process" << endl;
            // clearScreen();

            cin >> p2;
            if (p2 == 1)
            {
                Bisection();
                cout << endl;
                return 0;

                // 3
                // 1 0 -2 -5
            }
            else if (p2 == 2)
            {
                FalsePosition();
                cout << endl;
                return 0;

                // 3
                // 1 0 -2 -5
            }
            else if (p2 == 3)
            {
                Newtonraphson();
                cout << endl;
                return 0;

                // 3
                // 1 -6 11 -6
            }
            else if (p2 == 4)
            {
                Secant();
                cout << endl;
                return 0;

                // 3
                // 1 -6 11 -6
            }
            else if (p2 == 5)
            {
                cout << "The Process Is Terminated" << endl;
                return 0;
            }
            else
            {
                cout << "Invalid Input" << endl;
                return 0;
            }
        }
        else if (p == 3)
        {
            RungeKutta();
            return 0;
        }
        else if (p == 4)
        {
            NewtonForwardBackwardInterpolation();
            return 0;
        }
        else if (p == 5)
        {
            cout << "The Process Is Terminated" << endl;
            return 0;
        }
        else
        {
            cout << "Invalid Input" << endl;
            return 0;
        }
    }
}
