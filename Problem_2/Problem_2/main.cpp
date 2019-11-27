//  Problem_2
//  Created by Haoran Wang on 11/27/19.
//

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void getSpecialProducts (vector<int long> &v, int n) {
    
    int long left_prod[n], right_prod[n];
    
    left_prod[0] = 1;
    for (int i=1; i<n; i++)
        left_prod[i] = left_prod[i-1] * v[i-1];

    
    right_prod[n-1] = 1;
    for (int i=n-2; i>=0; i--)
        right_prod[i] = right_prod[i+1] * v[i+1];
    
    for (int i=0; i<n; i++)
        v[i] = left_prod[i] * right_prod[i];
    

}

int main(int argc, const char * argv[]) {

    ifstream input_file(argv[1]);
    
    if (!input_file)//check input files
    {
        cout << "Error: " << argv[1] << " could not be found!" << endl;
        exit(EXIT_FAILURE);
    }
    
    //Read input file
    int value;
    unsigned int n=0; //length of vector
    vector<int long> v;
    
    while (input_file >> value)
    {
        v.push_back(value);
        n++;
    }
    
    getSpecialProducts(v, n);
    
    //write output file
    ofstream myfile;
    myfile.open ("getSpecialProduct.txt");
    
    for (int i=0; i<n; i++)
    {
        myfile << v[i] << endl;
    }
    
    cout << "Job Completed! Please find results in getSpecialProduct.txt" << endl;
    
    return 0;
}
