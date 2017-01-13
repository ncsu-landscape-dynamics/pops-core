/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015 
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */


#include "Img.h"

using namespace std;


Img::Img(){	
	width = 0;
	height = 0;
	w_e_res = 0;
	n_s_res = 0;
	data  = NULL;
}

/*
Img::Img(int width,int height){
	this->width = width;
	this->height = height;
	data = (int *)std::malloc(sizeof(int)*width*height);
}
*/

Img::Img(int width, int height, int w_e_res,int n_s_res, int **data){
	this->width = width;
	this->height = height;
	this->w_e_res = w_e_res;
	this->n_s_res = n_s_res;
	this->data = data;
}

Img::Img(const char* fileName){
	GDALDataset *dataset;
	GDALRasterBand *dataBand;
	GDALAllRegister();
	dataset = (GDALDataset *)GDALOpen(fileName, GA_ReadOnly);
	double adfGeoTransform[6];
	if(!dataset) {
		cerr << "Can not open the image!" << endl;
	}else{
		width = dataset->GetRasterXSize();
		height = dataset->GetRasterYSize();

		if(dataset->GetGeoTransform(adfGeoTransform)==CE_None){
			w_e_res = abs(adfGeoTransform[1]);
			n_s_res = abs(adfGeoTransform[5]);
		}

		//cout << width << "x" << height <<endl;
		//cout << w_e_res << "X" << n_s_res << endl;

		dataBand = dataset->GetRasterBand(1);
		data = (int **)std::malloc(sizeof(int *) * height);
		int *stream = (int *)std::malloc(sizeof(int)*width*height);
		dataBand->RasterIO(GF_Read,0,0,width,height,stream,
			width,height,GDT_Int32,0,0);

		for(int i=0;i<height;i++){
			data[i] = &stream[i*width];
		}
		GDALClose((GDALDatasetH) dataset);
	}
}

int Img::getWidth(){
	return this->width;
}

int Img::getHeight(){
	return this->height;
}

int Img::getWEResolution(){
	return this->w_e_res;
}

int Img::getNSResolution(){
	return this->n_s_res;
}
/*
const int ** Img::getData(){
	return this->data;
}
*/


Img Img::operator+(Img& image){
	int re_width=0;
	int re_height=0;
	int **re_data=NULL;
	if(this->width != image.getWidth() || this->height != image.getHeight()){
		cerr << "The height or width of one image do not match with that of the other one!" << endl;
		return Img();
	}else{
		re_width = this->width;
		re_height = this->height;
		re_data = (int **)std::malloc(sizeof(int *) * re_height);
		int *stream = (int *)std::malloc(sizeof(int)*re_width*re_height);

		for(int i=0;i<re_height;i++){
			re_data[i] = &stream[i*re_width];
		}
		
		for(int i=0;i<re_height;i++){
			for(int j=0;j<re_width;j++){
				re_data[i][j] = this->data[i][j] + image.data[i][j];
			}
		}
		return Img(re_width,re_height,this->w_e_res, this->n_s_res,re_data);
	}
}

Img Img::operator-(Img& image){
	int re_width=0;
	int re_height=0;
	int **re_data=NULL;
	if(this->width != image.getWidth() || this->height != image.getHeight()){
		cerr << "The height or width of one image do not match with that of the other one!" << endl;
		return Img();
	}else{
		re_width = this->width;
		re_height = this->height;
		re_data = (int **)std::malloc(sizeof(int *) * re_height);
		int *stream = (int *)std::malloc(sizeof(int)*re_width*re_height);

		for(int i=0;i<re_height;i++){
			re_data[i] = &stream[i*re_width];
		}
		
		for(int i=0;i<re_height;i++){
			for(int j=0;j<re_width;j++){
				re_data[i][j] = this->data[i][j] - image.data[i][j];
			}
		}
		return Img(re_width,re_height,this->w_e_res, this->n_s_res,re_data);
	}
}

Img Img::operator*(int factor){

	int re_width = this->width;
	int re_height = this->height;
	int **re_data = (int **)std::malloc(sizeof(int *) * re_height);
	int *stream = (int *)std::malloc(sizeof(int)*re_width*re_height);

	for(int i=0;i<re_height;i++){
		re_data[i] = &stream[i*re_width];
	}
		
	for(int i=0;i<re_height;i++){
		for(int j=0;j<re_width;j++){
			re_data[i][j] = this->data[i][j] * factor;
		}
	}
	return Img(re_width,re_height,this->w_e_res, this->n_s_res,re_data);
}

Img::~Img(){
	if(data){
		if(data[0])
			delete [] data[0];
		delete [] data;
	}
}

