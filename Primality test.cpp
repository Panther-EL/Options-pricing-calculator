int number,count=0;
             cout <<"Enter your number:";
             cin>>number;
             for(int a=1; a<=number; a++){
                 if(number%a==0){
                     count++;
                 }
             }
             if (count==2){
                 cout<<"Number is prime";
             }
             else{
     cout<<"Number isn't prime";
       }
            
