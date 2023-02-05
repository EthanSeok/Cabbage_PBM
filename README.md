## Cabbage_PBM 2023

### Directory
* colab: google colab에서 실행할 수 있습니다. 해당 파일을 다운로드 한 후 구글 드라이브에 업로드 후 사용합니다.
* excel: excel을 이용하여 파라미터 관계 구현
* ipynb: Jupyter Notebook을 이용하여 각 함수 분리
* model: python cabbage model

### 배추 모형 소개
* https://docs.google.com/presentation/d/1mzTl4Q0dJcN_K1I_d30VK1Oqns8d3dC_/edit


### 강의영상
* https://youtube.com/playlist?list=PLqlg2QVJKrc244GkJrCjb0CoslLTns841


### 모델 분석 - 진행중
* https://docs.google.com/presentation/d/1VcwzRiALw2dmEhqGf25OCyT3Atbh8yvgGs2oLIkWt4A/edit


### 레퍼런스
* https://drive.google.com/drive/u/0/folders/1opr1q38xW6uWc002nPbcKC6waExuqCJo


### 배추 모형 파라미터 확인 - google colab
* https://drive.google.com/drive/u/0/folders/1nys0_ah4RvwPSNoYy0SheOZUPtwiWxba


<br>



## 모델 오류 수정사항
* fraction.py의 sun inclination 오류 정정


```
incl = np.arccos(sin_a + cos_b * np.cos(ha))
```

incl (sun inclination) 변수를 아래와 같이 수정함. 


```
incl = np.arcsin(sin_a + cos_b * np.cos(ha))
```


[de Fury, 1997](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-3040.1997.00094.x)에서 


<br>


다음 부분을 참고하였음.


<img src='https://user-images.githubusercontent.com/93086581/215049120-247c38fd-6a17-4e0c-a981-b64a841b9796.png'>


<img src='https://user-images.githubusercontent.com/93086581/215054014-36ed868f-1740-44b5-89e4-f93fa7516775.png'>


해당 예시를 참고하여 검증한 결과


**arccos** 사용시 
```
0.4863155364642975
```

**arcsin** 사용시
```
0.8737832677463231
```


**해당 결과를 통해 arcsin이 맞음을 검증.**

<br>


**arccos을 한 경우** 


<img src='https://user-images.githubusercontent.com/93086581/215050529-1d98c1d2-fb43-4a2c-9f5e-cad83b30cc35.png'>


**arcsin을 한 경우**


<img src='https://user-images.githubusercontent.com/93086581/215050700-05b4d99b-8459-4968-bd5d-d6b19ccf8616.png'>
