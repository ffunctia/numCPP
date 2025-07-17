#pragma once

#include <iostream>

using std::cout;
using std::endl;

namespace error{
    void errorCode(int errorNum){
        switch(errorNum){
            case 1:
            cout<<"\033[31m";
            cout<<"only the inverse of a square matrix can be found"<<"\033[0m";
            break;

            case 2:
            cout<<"\033[31m";
            cout<<"The inverse of a matrix with determinant 0 cannot be found."<<"\033[0m";
            break;

            case 3:
            cout<<"\033[31m";
            cout<<"Select a square matrix for adj()"<<"\033[0m";
            break;

            case 4:
            cout<<"\033[31m";
            cout<<"Select a square matrix for co_factor_matrix"<<"\033[0m";
            break;

            case 5:
            cout<<"\033[31m";
            cout<<"only the co_factor of a square matrix can be found"<<"\033[0m";
            break;

            case 6:
            cout<<"\033[31m";
            cout<<"only the minor of a square matrix can be found"<<"\033[0m";
            break;

            case 7:
            cout<<"\033[31m";
            cout<<"only the determinant of a square matrix can be found"<<"\033[0m";
            break;

            case 8:
            cout<<"\033[31m";
            cout<<"Dimensions for set_col do not match"<<"\033[0m";
            break;

            case 9:
            cout<<"\033[31m";
            cout<<"Dimensions for set_row do not match"<<"\033[0m";
            break;

            case 10:
            cout<<"\033[31m";
            cout<<"There is no such column index"<<"\033[0m";
            break;

            case 11:
            cout<<"\033[31m";
            cout<<"There is no such row index"<<"\033[0m";
            break;

            case 12:
            cout<<"\033[31m";
            cout<<"Dimensions for insert_column do not match"<<"\033[0m";
            break;

            case 13:
            cout<<"\033[31m";
            cout<<"Dimensions for insert_row do not match"<<"\033[0m";
            break;
            
            case 14:
            cout<<"\033[31m";
            cout<<"The number of columns of the first matrix must be equal to the number of rows of the other"<<"\033[0m";
            break;

            case 15:
            cout<<"\033[31m";
            cout<<"For orthogonal the dimensions must be equal"<<"\033[0m";
            break;

            case 16:
            cout<<"\033[31m";
            cout<<"For distance the dimensions must be equal"<<"\033[0m";
            break;

            case 17:
            cout<<"\033[31m";
            cout<<"For parallel the dimensions must be equal"<<"\033[0m";
            break;

            case 18:
            cout<<"\033[31m";
            cout<<"There is 0 in the vector, the application cannot be executed."<<"\033[0m";
            break;

            case 19:
            cout<<"\033[31m";
            cout<<"There is 0 in the vector, the application cannot be executed."<<"\033[0m";
            break;

            default:
            cout<<"error message undefined";

    };

    cout<<" ("<<errorNum<<")"<<endl;

    }
}