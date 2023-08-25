#include <iostream>

int r_binsearch(double val, double* arr, int len, int idx) {
    if (len==1) return idx;
    int mid = len/2;
    // arr has indices i i+1 ... i+len-1
    if (*(arr + mid) > val) {
        return r_binsearch(val, arr, len/2, idx);
    } else {
        return r_binsearch(val, arr+len/2, len - (len/2), idx+len/2);
    }
}
using namespace std;
int main() {
    double arr[] = {1,2,5,7,8,9,11,18,30};
    int i =r_binsearch(9.1, arr, 9, 0);
    cout<<i<<" "<<arr[i]<<endl;
    i =r_binsearch(10.9, arr, 9, 0);
    cout<<i<<" "<<arr[i]<<endl;
    i =r_binsearch(-1, arr, 9, 0);
    cout<<i<<" "<<arr[i]<<endl;
    i =r_binsearch(11.1, arr, 9, 0);
    cout<<i<<" "<<arr[i]<<endl;
    i =r_binsearch(1.1, arr, 9, 0);
    cout<<i<<" "<<arr[i]<<endl;
}
