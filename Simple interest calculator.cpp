 //Assignment 5(Simple Intrest Calculator)
#include <iostream>
using namespace std;
int main(){
     int time;
     cout << "Enter the time:";
     cin >>time;
     int principal;
     cout << "Enter the principal:";
     cin >>principal;
     int rate;
     cout << "Enter the rate:";
     cin >> rate;
     int interest =(principal*rate*time)/100;
     cout << interest;
}
