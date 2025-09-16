#include <iostream>
using namespace std;
int main(){
int number;
int sum = 0;
cout <<"Enter your number:"<<endl;
cin >>number;
for (int i = 1; i<number; i++){
sum=sum+i;
}
cout<<"The sum from 1 to"<<number<<"is:"<<sum<<endl;
}
