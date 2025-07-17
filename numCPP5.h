#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "errors.h"
#include <iomanip>

using std::vector;
using std::cout;
using std::endl;
using error::errorCode;
using std::unique_ptr;
using std::make_unique;
using std::move;
using std::swap;
using std::setprecision;
using std::setw;
using std::fixed;
using std::abs;

namespace hidden{
    bool isZeroElement(vector<double> vecForisZero){
        for(int i=0; i<vecForisZero.size(); i++){
            if(abs(vecForisZero[i]) < 1e-10){
                return true;
            }
        }

        return false;
    }
}

using namespace hidden;

namespace vec_tools{

    void print_vector(const vector<double>& vecForPrint){
        for(int i=0; i<vecForPrint.size(); i++){
            cout<<vecForPrint[i]<<", ";
        }
        cout<<endl;
    }

    double magnitute(const vector<double>& vecForMagnitute){
        double vecMagnitute =0;
        for(int i = 0; i< vecForMagnitute.size(); i++){
            vecMagnitute+=vecForMagnitute[i] * vecForMagnitute[i];
        }
        vecMagnitute = sqrt(vecMagnitute);

        return vecMagnitute;
    }

    double scalar_dot(const vector<double>& vecForScalarDot1,
    const vector<double>& vecForScalarDot2){
        if(vecForScalarDot1.size() != vecForScalarDot2.size()){
            return NAN;
        }
        double scalarDot = 0;
        for(int i = 0; i<vecForScalarDot1.size(); i++){
            scalarDot+= vecForScalarDot1[i] * vecForScalarDot2[i];
        }

        return scalarDot;
    }

    double angle(const vector<double>& vecForAngle1,
    const vector<double>& vecForAngle2){
        double magnitute1 = magnitute(vecForAngle1);
        double magnitute2 = magnitute(vecForAngle2);
        if(isZeroElement(vecForAngle1) || isZeroElement(vecForAngle2)){
            errorCode(19);
            return NAN;
        }
        double scalarDotofVectors = scalar_dot(vecForAngle1,vecForAngle2);

        double angle = acos(scalarDotofVectors/(magnitute1*magnitute2));

        return angle;

    }

    double angle_deg(const vector<double>& vecForAngle1, 
        const vector<double>& vecForAngle2){
            return angle(vecForAngle1, vecForAngle2) * 180.0 / M_PI;
    }


    vector<double> cross_product(const vector<double>& vecForDot1,
    const vector<double>& vecForDot2){
        if(vecForDot1.size() != vecForDot2.size() || vecForDot1.size()!=3){
            return vecForDot1;
        }

        vector<double> dottedVec;
        dottedVec.push_back(vecForDot1[1]*vecForDot2[2] - vecForDot1[2]*vecForDot2[1]);
        dottedVec.push_back(vecForDot1[2]*vecForDot2[0] - vecForDot1[0]*vecForDot2[2]);
        dottedVec.push_back(vecForDot1[0]*vecForDot2[1] - vecForDot1[1]*vecForDot2[0]);

        return dottedVec;

    }

    vector<double> element_wise_dot(const vector<double>& vecForDot1, 
        const vector<double>& vecForDot2){
            if(vecForDot1.size() != vecForDot2.size()) return vecForDot1; 
            vector<double> result;
            for(int i = 0; i < vecForDot1.size(); i++){
                result.push_back(vecForDot1[i]*vecForDot2[i]);
            }
            return result;
    }

    vector<double> vector_sum(const vector<double>& vecForSum1,
    const vector<double>& vecForSum2){
        if(vecForSum1.size() != vecForSum2.size()){
            return vecForSum1;
        }

        vector<double> sumVector;
        for(int i=0; i<vecForSum1.size(); i++){
            sumVector.push_back(vecForSum1[i]+vecForSum2[i]);
        }

        return sumVector;
    }

    vector<double> vector_sub(const vector<double>& vecForSub1,
    const vector<double>& vecForSub2){
        if(vecForSub1.size() != vecForSub2.size()){
            return vecForSub1;
        }

        vector<double> subVector;
        for(int i=0; i<vecForSub1.size(); i++){
            subVector.push_back(vecForSub1[i]-vecForSub2[i]);
        }

        return subVector;
    }

    vector<double> multi(const vector<double>& vecForMulti, double num){
        vector<double> returnedVector;
        for(int i=0; i<vecForMulti.size(); i++){
            returnedVector.push_back(vecForMulti[i]*num);
        }

        return returnedVector;
    }

    vector<double> convert_unit(const vector<double>& vecForUnit){
        vector<double> returnedVector;
        double norm = magnitute(vecForUnit);
        if(norm==0){
            errorCode(1);
            return vecForUnit;
        }
        returnedVector = multi(vecForUnit, 1/norm);

        return returnedVector;
    }

    vector<double> projection(const vector<double>& vecForPro1, 
        const vector<double>& vecForPro2){
            double scalar = scalar_dot(vecForPro1, convert_unit(vecForPro2));
            return multi(convert_unit(vecForPro2), scalar);
    }

    bool parallel(const vector<double>& vecForPar1,
    const vector<double>& vecForPar2){
        if(vecForPar1.size() != vecForPar2.size()){
            errorCode(17);
            return NAN;
        }

        if(isZeroElement(vecForPar1) || isZeroElement(vecForPar2)){
            errorCode(18);
            return NAN;
        }

        bool isParallel = true;
        double k=vecForPar1[0]/vecForPar2[0];
        for(int i=1; i<vecForPar1.size(); i++){
            if(vecForPar1[i]/vecForPar2[i] != k){
                isParallel = false;
                break;
            }
        }

        return isParallel;
    }

    double distance_between(const vector<double>& vecForDistance, 
        const vector<double>& vecForDistance2){
            if(vecForDistance.size() != vecForDistance2.size()){
                errorCode(16);
                return NAN;
            }
            double distance = 0.0;
            distance = magnitute(vector_sub(vecForDistance, vecForDistance2));
            return distance;
    }

    vector<double> normalize(const vector<double>& vecForNorm, double len){
        return multi(convert_unit(vecForNorm), len);
    }

    vector<double> zeros(int len){
        vector<double> zeroVec;
        for(int i=0; i<len; i++){
            zeroVec.push_back(0.0);
        }

        return zeroVec;
    }

    bool orthogonal(const vector<double>& vecForO1, const vector<double>& vecForO2){
        if(vecForO1.size() != vecForO2.size()){
            errorCode(15);
            return NAN;
        }
        else if(abs(scalar_dot(vecForO1, vecForO2)) < 1e-10){
            return true;
        }
        else{
            return false;
        }
    }

    bool zero_vector(const vector<double>& vecForZero){
        for(auto x:vecForZero){
            if(x!=0.0) return false;
        }

        return true;
    }  

}

using namespace vec_tools;
//classes
class matrix{
    public:
    unique_ptr<vector<vector<double>>> _matrix;
    unique_ptr<int> _nRow;
    unique_ptr<int> _nCol;
    unique_ptr<int> _sign;

    public:
    //Constructors
    matrix(){
        _matrix = make_unique<vector<vector<double>>>();
        _nRow = make_unique<int>();
        _nCol = make_unique<int>();
        _sign = make_unique<int>(1);
    }
    
    matrix(vector<double> vecForC, int nRow, int nCol){
        _matrix = make_unique<vector<vector<double>>>();

        for(int i = 0; i<nRow; i++){
            vector<double> row;
            for(int j=i; j<nRow*nCol; j+=nRow){
                row.push_back(vecForC[j]);
            }
            _matrix->push_back(row);
        }

        _nCol = make_unique<int>(nCol);
        _nRow = make_unique<int>(nRow);
        _sign = make_unique<int>(1);

    }

    matrix(const vector<vector<double>>& vecForC){
        _nRow = make_unique<int>(vecForC.size());
        _nCol = make_unique<int>(0);
        _sign = make_unique<int>(1);

        vector<vector<double>> vecForCC= vecForC;
        for(int i = 0; i<*_nRow; i++){
            if(vecForCC[i].size() > *_nCol){
                *_nCol = vecForCC[i].size();
            }
        }

        for(int i=0; i<*_nRow; i++){
            for(int j=vecForCC[i].size(); j<*_nCol; j++){
                vecForCC[i].push_back(0);
            }
        }

        _matrix = make_unique<vector<vector<double>>>(vecForCC);
    }

    matrix(const matrix& other){
        _matrix = make_unique<vector<vector<double>>>(*other._matrix);
        _nCol = make_unique<int>(*other._nCol);
        _nRow = make_unique<int>(*other._nRow);
        _sign = make_unique<int>(*other._sign);
    }

    matrix& operator=(const matrix& other){
        if (this != &other) {
            _matrix = make_unique<vector<vector<double>>>(*other._matrix);
            _nCol = make_unique<int>(*other._nCol);
            _nRow = make_unique<int>(*other._nRow);
            _sign = make_unique<int>(*other._sign);
        }
        return *this;
    }

    matrix(matrix&& other) noexcept{
        _matrix = move(other._matrix);
        _sign = move(other._sign);
        _nRow = move(other._nRow);
        _nCol = move(other._nCol);
    }

    matrix& operator=(matrix&& other) noexcept{
        if (this != &other) {
            _matrix = move(other._matrix);
            _sign = move(other._sign);
            _nRow = move(other._nRow);
            _nCol = move(other._nCol);
        }
        return *this;
    }

    //Functions

    void print_matrix(){
        for(int i = 0; i < *_nRow; i++){
            cout<<"|";
            for(int j = 0; j < *_nCol; j++){
                cout << fixed << setprecision(5) << setw(13) << (*_matrix)[i][j] << " ";
            }
            cout <<"|"<<endl;
        }
    }

    matrix transpose() const{
        
        int nCol = *_nRow;
        int nRow = *_nCol;
        
        vector<double> forTransposed;
        
        for(int i=0; i<*_nRow; i++){
            for(int j=0; j<*_nCol; j++){
                forTransposed.push_back((*_matrix)[i][j]);
            }
        }

        matrix transposedMatrix(forTransposed,nRow,nCol);

        return transposedMatrix;

    }

    matrix dot(const matrix& otherMatrix){
        if(*_nCol != *otherMatrix._nRow){
            errorCode(14);
            return *this;
        }

        matrix dotMatrix;
        matrix matrixForKeep = otherMatrix.transpose();

        for(int i=0; i<*_nRow; i++){

            vector<double> row;
            for(int j=0; j<*matrixForKeep._nRow; j++){
                double value = 0;
                
                for(int k=0; k<*_nCol; k++){
                    value += (*_matrix)[i][k] * (*matrixForKeep._matrix)[j][k];
                }
                
                row.push_back(value);
            }
        
            dotMatrix._matrix->push_back(row);

        }
        *dotMatrix._nRow = *_nRow;
        *dotMatrix._nCol = *otherMatrix._nCol;

        return dotMatrix;
    }

    matrix dot(double num){
        matrix dotMatrix(*this);

        for(int i=0; i<*_nRow; i++){
            for(int j=0; j<*_nCol; j++){
                (*dotMatrix._matrix)[i][j] *=num;
            }
        }

        return dotMatrix;

    } 

    matrix power(int times){
        matrix matrixForPower(*this);

        if(times == 0 ){
            return *this;
        }

        if(times <0){
            return this->inverse().power(-times);
        }

        for(int i = 0; i<times; i++){
            matrixForPower = matrixForPower.dot(*this);
        }

        return matrixForPower;
    }

    vector<double> get_row(int nRow) const{
        if(nRow >= *_nRow || nRow <0){
            errorCode(11);
            return {-999};
        }

        return (*_matrix)[nRow];
    }

    vector<double> get_col(int nCol) const{
        if(nCol >= *_nCol || nCol <0){
            errorCode(10);
            return {-999};
        }

        vector<double> col;
        for(int i=0; i<*_nRow; i++){
            col.push_back((*_matrix)[i][nCol]);
        }

        return col;
    }

    void insert_row(const vector<double>& row, int indexR){
        if(indexR >= *_nRow || indexR <0){
            errorCode(11);
            return;
        }

        if(row.size() != *_nCol){
            errorCode(13);
            return;
        }

        vector<vector<double>> forRowIns;
        forRowIns.push_back(row);
        _matrix->insert(_matrix->begin()+indexR, forRowIns.begin(), forRowIns.end());
        *_nRow+=1;
    }

    void insert_col(const vector<double>& col, int indexC){
        if(indexC >= *_nCol || indexC <0){
            errorCode(10);
            return;
        }

        if(col.size() != *_nRow){
            errorCode(12);
            return;
        }

        for(int i=0; i<*_nRow; i++){
            (*_matrix)[i].insert((*_matrix)[i].begin()+indexC, col[i]);
        }

        *_nCol+=1;
    }

    void delete_row(int indexR){
        if(indexR >= *_nRow || indexR <0){
            errorCode(10);
            return;
        }
        _matrix->erase(_matrix->begin()+indexR);
        *_nRow-=1;
    }

    void delete_col(int index){
        if(index >= *_nCol || index <0){
            errorCode(10);
            return;
        }
        for(int i=0; i<*_nRow; i++){
            (*_matrix)[i].erase((*_matrix)[i].begin()+index); 
        }

        *_nCol-=1;
    }

    void set_row(vector<double> vecForSet, int indexR){
        if(vecForSet.size() != *_nCol){
            errorCode(9);
            return;
        }

        (*_matrix)[indexR] = vecForSet;
    }

    void set_col(vector<double> vecForSet, int indexC){
        if(vecForSet.size() != *_nRow){
            errorCode(8);
            return;
        }
        for(int i =0; i<*_nRow; i++){
            (*_matrix)[i][indexC] = vecForSet[i];
        }
    }

    void swap_row(int index, int to){
        vector<double> memory = get_row(to);
        set_row(get_row(index), to);
        set_row(memory, index);
        *_sign *= -1;
    }

    void swap_col(int index, int to){
        vector<double> memory = get_col(to);
        set_col(get_col(index), to);
        set_col(memory, index);
        *_sign *= -1;
    }

    vector<int> zero_row_index(){
        vector<int> zeroRows;
        for(int i=0; i<*_nRow; i++){
            if(zero_vector(get_row(i))){
                zeroRows.push_back(i);
            }
        }

        return zeroRows;
    } 
    
    vector<int> zero_col_index(){
        vector<int> zeroCols;
        for(int i=0; i<*_nCol; i++){
            if(zero_vector(get_col(i))){
                zeroCols.push_back(i);
            }
        }

        return zeroCols;
    }

    double det() const{
        if(*_nCol != *_nRow){
            errorCode(7);
            return -999;
        }

        double determinant=1;
        matrix gauss = gauss_elimination();
        for(int i=0; i<*gauss._nRow; i++){
            if(abs((*gauss._matrix)[i][i]) < 1e-10){
                return 0;
            }
            determinant*=(*gauss._matrix)[i][i];
        }
        return *gauss._sign * determinant;
    }

    matrix minor_matrix(int nRow, int nCol) const{
        matrix returnedMatrix(*this);
        returnedMatrix.delete_col(nCol);
        returnedMatrix.delete_row(nRow);
        return returnedMatrix;
    }

    double minor(int nRow, int nCol) const{
        if(*_nCol != *_nRow){
            errorCode(6);
            return -999;
        }
        
        return minor_matrix(nRow, nCol).det();
    }

    double co_factor(int nRow, int nCol) const{
        if(*_nCol != *_nRow){
            errorCode(5);
            return -999;
        }
        int sign=-1;
        if((nRow+nCol)%2==0) sign=1;
        return minor(nRow, nCol) * sign;
    }

    matrix co_factor_matrix() const{
        if(*_nCol != *_nRow){
            errorCode(4);
            return *this;
        }
        return adj().transpose();
    }

    matrix adj() const{
        if(*_nCol != *_nRow){
            errorCode(3);
            return *this;
        }
        return inverse().dot(det());
    }

    matrix inverse() const{
        if(*_nCol != *_nRow){
            errorCode(1);
            return *this;
        }

        if(det() == 0){
            errorCode(2);
            return *this;
        }

        matrix matrixForReturn(*this);
        for(int i=0; i<*_nRow; i++){
            for(int j=0; j<*_nCol; j++){
                if(j==i){
                    (*matrixForReturn._matrix)[i].push_back(1);
                    continue;
                }

                (*matrixForReturn._matrix)[i].push_back(0);
            }
        }

        *matrixForReturn._nCol *=2;
        matrixForReturn = matrixForReturn.gauss_jordan_elimination();

        for(int i=0; i<*_nCol; i++){
            matrixForReturn.delete_col(0);
        }

        return matrixForReturn;
        
    }

    matrix gauss_elimination() const{
        matrix returnedMatrix(*this);
        vector<int> zeroCol = returnedMatrix.zero_col_index();

        for(int i=0; i<*returnedMatrix._nRow; i++){
            bool isZero=false;
            for(int j=0; j<zeroCol.size(); j++){
                if(i==zeroCol[j]){
                    isZero=true;
                    break;
                }
            }

            if(isZero) continue;

            if(abs((*returnedMatrix._matrix)[i][i]) <1e-10){
                for(int j=i+1; j<*returnedMatrix._nRow; j++){
                    if((*returnedMatrix._matrix)[j][i]!=0){
                        returnedMatrix.swap_row(i,j);
                        break;
                    }
                }
            }

            if(abs((*returnedMatrix._matrix)[i][i]) <1e-10) continue;

            for(int j = i+1; j<*returnedMatrix._nRow; j++){
                double factor = (*returnedMatrix._matrix)[j][i] / (*returnedMatrix._matrix)[i][i];
                (*returnedMatrix._matrix)[j] = vector_sub((*returnedMatrix._matrix)[j],
                multi((*returnedMatrix._matrix)[i],factor));

            }
        }

        return returnedMatrix;
        
    }

    matrix gauss_jordan_elimination() const {
        matrix returnedMatrix(*this);

        int lead = 0;
        for (int r = 0; r < *returnedMatrix._nRow; ++r) {
            if (lead >= *returnedMatrix._nCol)
                break;

            int i = r;
            while (i < *returnedMatrix._nRow && abs((*returnedMatrix._matrix)[i][lead]) < 1e-10)
                ++i;

            if (i == *returnedMatrix._nRow) {
                ++lead;
                --r; 
                continue;
            }

            swap((*returnedMatrix._matrix)[i], (*returnedMatrix._matrix)[r]);

            double pivot = (*returnedMatrix._matrix)[r][lead];
            for (int j = 0; j < *returnedMatrix._nCol; ++j)
                (*returnedMatrix._matrix)[r][j] /= pivot;

            for (int i = 0; i < *returnedMatrix._nRow; ++i) {
                if (i != r) {
                    double factor = (*returnedMatrix._matrix)[i][lead];
                    for (int j = 0; j < *returnedMatrix._nCol; ++j)
                        (*returnedMatrix._matrix)[i][j] -= factor * (*returnedMatrix._matrix)[r][j];
                }
            }

            ++lead;
        }

        return returnedMatrix;
    }

    int rank() const{
        matrix matrixForRank = gauss_jordan_elimination();
        int rankValue = 0;
        for(int i=0; i<*_nRow; i++){
            vector<double> row = matrixForRank.get_row(i);
            if(zero_vector(row)) continue;
            else rankValue++;
        }

        return rankValue;
    }
            
};

