# Dual-Band FMCW Radar Fusion Algorithms (MATLAB Implementation)

This repository contains MATLAB implementations of **dual-band radar fusion algorithms** developed for a low-cost FMCW radar system.  
The fusion is performed **at the signal level** (after dechirping, before FFT), enabling improved target detection under low-SNR conditions.

<p align="center">
  <img src="pics/FMCW E.drawio.png"/>
</p>

## 📌 Implemented Fusion Algorithms
- **GEF:** Gain Envelope Fusion  
- **AEF:** Adaptive Envelope Fusion  
- **GPF:** Gain Power Fusion  
- **APF:** Adaptive Power Fusion  

---

## 📚 Theoretical Background

### 🔹 Fusion Stage
The fusion is applied on the **beat signals** `s1(t)` and `s2(t)` from the 5.8 GHz and 24 GHz radars,  
after bandpass filtering and envelope/power extraction.

---

### 🔹 Envelope Fusion (GEF)
<p align="center">
  <img src="pics/eq1.svg" width="300"/>
</p>

<p align="left"><i>where w1 and w2 are fixed fusion weights.</i></p>

---

### 🔹 Adaptive Envelope Fusion (AEF)
<p align="center">
  <img src="pics/eq2.svg" width="320"/>
</p>

<p align="center">
  <img src="pics/eq3.svg" width="250"/>
</p>

<p align="left"><i>where Q_i(t) is the segment energy of radar i.</i></p>

---

### 🔹 Power Fusion (GPF)
<p align="center">
  <img src="pics/eq4.svg" width="300"/>
</p>

---

### 🔹 Adaptive Power Fusion (APF)
<p align="center">
  <img src="pics/eq5.svg" width="320"/>
</p>

<p align="center">
  <img src="pics/eq6.svg" width="250"/>
</p>

<p align="left"><i>where SNR_i(t) is the instantaneous signal-to-noise ratio of radar i.</i></p>

---

## 📂 Repository Structure
```
├── pics/                  # Equation images for README
│ ├── eq1.svg
│ ├── eq2.svg
│ ├── eq3.svg
│ ├── eq4.svg
│ ├── eq5.svg
│ └── eq6.svg
├── data/                  # Example radar beat signal files
│ ├── beat_5_8GHz.mat
│ └── beat_24GHz.mat
├── GEF_Fusion.m           # MATLAB script for Global Envelope Fusion
├── AEF_Fusion.m           # MATLAB script for Adaptive Envelope Fusion
├── GPF_Fusion.m           # MATLAB script for Global Power Fusion
├── APF_Fusion.m           # MATLAB script for Adaptive Power Fusion
└── README.md              # Project documentation
```
---

## 🚀 How to Use

### ✅ Requirements
- MATLAB R2020a or later
- Signal Processing Toolbox

### ✅ Steps to Run
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/Radar-Fusion.git
   cd Signal-Fusion-Dual-Band-FMCW-Radar
   ```
2. Open MATLAB and navigate to the project folder.
3. Load the provided example beat signals (data/beat_5_8GHz.mat, data/beat_24GHz.mat).
4. Run any of the fusion scripts:
   ```
   Run any of the fusion scripts:
   ```

The script will:

1. Loads radar beat signals.
2. Normalizes and processes the signals.
3. Applies the selected fusion algorithm.
4. Estimates the fused beat frequency and range.
5. Plots the input and fused signals.

## 📊 Example Output
The scripts generate:
1. Time-domain plots of the individual and fused beat signals.
2. FFT spectra of the fused signals.
3. Plots of adaptive weights for AEF and APF.
```
fb_5.8GHz   = 37.354 MHz
fb_24GHz    = 37.500 MHz
fb_fused    = 37.420 MHz
R_fused     = 74.8 m
```

## 📝 Citation
If you use this code in your research, please cite:
@article{YourName2025RadarFusion,
  title={Non-Coherent Power and Envelope Fusion Algorithms for Low-Cost Dual-Band FMCW Radar},
  author={Your Name},
  journal={To be submitted},
  year={2025}
}

## 📄 License
This project is licensed under the MIT License – you are free to use, modify, and distribute it with proper attribution.
---
