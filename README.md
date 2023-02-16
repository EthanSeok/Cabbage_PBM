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

---


### 모델 오류 수정사항
#### fraction.py의 sun inclination 수정사항


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
`0.4863155364642975`


**arcsin** 사용시
`0.8737832677463231`


**해당 결과를 통해 arcsin이 보다 적합하다고 판단.**

<br>


**(좌) arccos을 한 경우 (24시간)**    **(우) arcsin을 한 경우 (24시간)**

![image.jpg1](https://user-images.githubusercontent.com/93086581/218964817-eae2d83d-eecf-4263-baf1-90629222d981.png) |![image.jpg2](https://user-images.githubusercontent.com/93086581/218964827-b5539969-6ba5-4998-8ccf-30ed955fd269.png)
--- | --- | 



<br>

**하지만 아래 그래프와 같이 arccos를 arcsin으로 변경시 sunhgt 값이 적합하지 않은 결과를 보임.**


<img src='https://user-images.githubusercontent.com/93086581/218955037-11390f35-b07d-4811-bfce-07cb0a401775.png'>
 
 
<br>


확인한 결과 arccos을 사용했을 때의 sunhgt와 arcsin을 사용했을 때의 incl이 동일함을 확인.

**(좌) arccos 사용시 sunhgt**, **(중) arcsin 사용시 incl**, **(우) arcsin 사용 및 max(0.05) 사용시 incl**
![image.jpg1](https://user-images.githubusercontent.com/93086581/218956957-89c3ca4a-dfe1-4fd6-8222-90c3b68e3944.png) |![image.jpg2](https://user-images.githubusercontent.com/93086581/218956977-ba07b0e5-b767-44d7-af77-444587b413ea.png) |![image.jpg3](https://user-images.githubusercontent.com/93086581/218957684-0f3faeb8-3eb4-4aa3-bf69-98ba47f5638a.png)
--- | --- | --- |

<img src='https://user-images.githubusercontent.com/93086581/218959176-b1d1eea3-5086-4647-a541-31dfde742718.png'>

따라서 최종적으로 다음과 같이 수정함.


```
sunhgt = np.arcsin(max(0.05(sin_a + cos_b * np.cos(ha)))
```

---


#### stage.py의 calcVerdvs 수정사항

```
    def calcVerdvs(self, Ta):
        Ta  = max(Ta, 0.01)
        rate = np.exp(-1*(np.log(Ta/optVer)**4))
        self.sumVer += rate * conv
        self.verdvs = max(1, self.sumVer/satVer)
```


를 아래와 같이 수정


```
    def calcVerdvs(self, Ta):
        Ta  = max(Ta, 0.01)
        rate = np.exp(-1*(np.log(Ta/optVer)**4))
        self.sumVer += rate * conv
        self.verdvs = min(1, self.sumVer/satVer)
```


<img src='https://user-images.githubusercontent.com/93086581/218971875-cc48db2b-d6bf-4396-9fcc-44b10b1cbe42.png'>


* `self.verdvs = max(1, self.sumVer/satVer)` 부분에서 max를 사용할 시 verdvs의 최소가 1이 되므로 bolting 시기가 매우 앞당겨짐. 따라서 이를 min으로 수정하여 bolting 시기 정정


* max 사용시 DW 및 FW의 결과가 매우 낮았음. min으로 수정하니, 보다 합리적인 결과가 도출됨.

* 해당 결과를 통해 min 사용이 더 적합하다고 판단.

### 챔버환경 적용사항
* 챔버환경 적용을 위한 fraction.py의 광 환경을 다음과 같이 조정함.

```
sunhgt = np.arcsin(max(0.05,sin_a + cos_b * np.cos(ha)))
        if sunhgt > 0.051: sunhgt = 1
```

<img src = 'https://user-images.githubusercontent.com/93086581/219265396-5feb0016-f202-406a-88b6-eff55d4a004c.png'>

* 광량에 따른 LAI의 변화 ==> An 값 변화로 인하여 DW 및 FW에만 변가 있었음.
* 챔버 환경인 300µmol로 고정시 DW및 FW 저하
