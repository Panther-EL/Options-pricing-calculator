#include <iostream>
#include <string>
using namespace std;


class Student{
    private:
     int id;
     string name;
     float scores[5];

     public:
      void inputDetails();
      void displayDetails() const;
      float calculateAverage() const;
      char calculateGrade() const;
      int getId() const;

};

void Student::inputDetails(){
    cout << "Enter student ID:";
    cin >> id;
    cin.ignore();
    cout << "Enter student's name:";
    getline(cin,name);
    cout << "Enter the five(5) scores of the student:" <<endl;
    for(int i =0; i<=4; i++){
        cout <<"Score" << i+1 <<":";
        cin >> scores[i];
    }
}



void Student::displayDetails()const{
    cout <<"Student ID:"<<id<<endl;
    cout <<"Name:"<<name<<endl;
     cout << "Scores:";
     for(int i =0; i<=4; i++){
         cout << scores[i]<< " " <<endl;
     }
    cout <<"Average:" <<calculateAverage()<<endl;
    cout <<"Grade:" <<calculateGrade()<<endl;
}


float Student::calculateAverage()const{
    float sum = 0;
    for(int i=0; i<5; i++){
        sum += scores[i];
    }
    return sum/5;
}


char Student::calculateGrade()const{
    float avg = calculateAverage();
    if(avg>=80) return 'A';
    else if(avg>=70) return 'B';
    else if(avg>=60) return 'C';
    else if(avg>=50) return 'D';
    else return 'F';
}

int main(){
    Student s;
    s.inputDetails();
    s.displayDetails();
