void FalsePosition()
// {
//     int p;
//     cout << "Enter power: ";
//     cin >> p;

//     vector<int> arr(p + 1);
//     cout << "Enter coefficients: ";
//     for (int i = 0; i < arr.size(); i++)
//     {
//         cin >> arr[i];
//     }

//     // double xMax = ceil(sqrt(pow(arr[1] / (double)arr[0], 2) - 2 * (arr[2] / (double)arr[0])));
//     // double r1 = -xMax;
//     // double r2 = xMax;
//     double r1=-1000;
//     double r2=1000;

//     double a = r1, b = r2;
//     bool intervalFound = false;

//     for (int i = r1; i <= r2; i++)
//     {
//         if ((fx(i, arr) * fx(i + 1, arr)) < 0)
//         {
//             a = i;
//             b = i + 1;
//             intervalFound = true;
//             break;
//         }
//     }

//     // if (!intervalFound || fx(a, arr) * fx(b, arr) >= 0)
//     // {
//     //     cout << "Enter a correct range where the function changes sign." << endl;
//     //     return;
//     // }

//     double tolerance = 1e-6;
//     int max_iterations = 1000;
//     double c, fa = fx(a, arr), fb = fx(b, arr);

//     for (int i = 1; abs(b - a) >= tolerance && i <= max_iterations; i++)
//     {
//         c = (a * fb - b * fa) / (fb - fa);

//         if (fabs(fx(c, arr)) < tolerance)
//         {
//             cout << "Exact root is found within tolerance." << endl;
//             cout << "root = " << c << endl;
//             return;
//         }

//         double fc = fx(c, arr);

//         if (fa * fc < 0)
//         {
//             b = c;
//             fb = fc; // Update fb to fx(b, arr)
//         }
//         else
//         {
//             a = c;
//             fa = fc; // Update fa to fx(a, arr)
//         }
//     }

//     cout << "root = " << c << endl;
// }