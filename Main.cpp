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
                if (j < i)
                {
                    sum -= AugmentedMatrix[i][j] * x[j];
                }
                else
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


int main()
{
    int p;
    while (1)
    {
        cout << "Press 1 to solve system of linear equations" << endl;
        cout << "press 2 to terminate the process" << endl;
        cin >> p;
        clearScreen();

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
