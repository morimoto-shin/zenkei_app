import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

st.set_page_config(page_title='惑星会合シミュレーター')
st.markdown("<h1 style='text-align: center; color: black;'>惑星会合シミュレーター</h1>", unsafe_allow_html=True)
input_elev = st.sidebar.slider("xy軸からの方位角", 0, 180, 45)
input_azim =  st.sidebar.slider("z軸からの方位軸", 0, 180, 45)

#"""----------関数----------"""
# 楕円軌道の関数
def ellipse(i, a, e):
    u=np.pi/180*i
    b=a*np.sqrt(1-np.power(e,2))
    x=a*np.cos(u)-a*e
    y=b*np.sin(u)
    return x, y

# 楕円軌道の関数
def ellipse_Kepller(u, a, e):
    b=a*np.sqrt(1-np.power(e,2))
    x=a*np.cos(u)-a*e
    y=b*np.sin(u)
    return x, y

# 日本標準時（JST）からユリウス暦（JD）への変換関数
def DAYtoJD(Year,Month,Day,Hour,Minute):
    DDay = Day+(Hour+Minute/60-9)/24
    if Month < 3:
        MMonth = Month +12
        YYear = Year - 1
    else:
        MMonth = Month
        YYear = Year
    aa = int(YYear/100)
    bb = 2.0 - aa + int(aa/4.0)
    JD = int(365.25*YYear) + int(30.6001*(MMonth+1)) + DDay + bb + 1720994.5
    return JD

# ユリウス暦（JD）から日本標準時（JST）への変換関数
def JDtoDAY(JD):
    JDs = JD + 0.5 + 9/24
    JDi = int(JDs)
    JDf = JDs - JDi
    if JDi > 2299160.0:
        aa = int((JDi - 1867216.25) / 36524.25)
        JDi = JDi + (1.0 + aa - int(aa / 4.0))
    bb = JDi + 1524.0
    YY = int((bb - 122.1) / 365.25)
    DD = bb - int(YY * 365.25)
    MM = int(DD / 30.6001)
    DDay= DD - int(MM * 30.6001) + JDf
    Day = int(DDay)
    HH = DDay - Day
    Hour = int(HH*24)
    Minute = round((HH*24 - Hour)*60)
    if MM < 14:
        Month = MM - 1
    else:
        Month = MM - 13
    if Month > 2:
        Year = YY - 4716
    else:
        Year = YY - 4715
    return Year, Month, Day, Hour, Minute

# ケプラー方程式の関数
def funct(u,t,e,mm,t0):
    return u - e*np.sin(u) - mm*(t - t0)

# ケプラー方程式の関数の微分
def diff_funct(u,e):
    return 1 - e*np.cos(u)

# 座標変換（軌道面→日心直交座標）
def coordtrans(xx,yy,cosohm,sinohm,cosi,sini,cosomega,sinomega):
    Xc=xx*(cosohm*cosomega-sinohm*cosi*sinomega)-yy*(cosohm*sinomega+sinohm*cosi*cosomega)
    Yc=xx*(sinohm*cosomega+cosohm*cosi*sinomega)-yy*(sinohm*sinomega-cosohm*cosi*cosomega)
    Zc=xx*sini*sinomega+yy*sini*cosomega
    return Xc,Yc,Zc

def main():
    #"""----------惑星クラス----------"""
    class Planet():
        def __init__(self, a, e, ohm, i, omega, T, b, n, t, color):
            self.a = a  # 軌道長半径
            self.e = e  # 軌道離心率
            self.ohm = ohm #
            self.i = i
            self.omega = omega
            self.T = T
            self.b = b 
            self.n = n
            self.t = t
            self.cosohm = np.cos(ohm)
            self.sinohm = np.sin(ohm)
            self.cosi = np.cos(i)
            self.sini = np.sin(i)
            self.cosomega = np.cos(omega)
            self.sinomega = np.sin(omega)
            self.color = color

    #"""----------各惑星のインスタンス生成----------"""
    planets = []

    # 地球
    planets.append(Planet(
        1.00000102, # a
        0.01671022, # e
        np.pi/180* 174.838, # ohm
        np.pi/180* 0.002, # i
        np.pi/180*(102.972 - 174.838), # omega
        365.24, # T
        1.00000102*np.sqrt(1-np.power(0.01671022,2)), # b
        2*np.pi/365.24, # n
        2455199.50625, # t
        'blue',
    ))

    # 火星
    planets.append(Planet(
        1.52371034,
        0.09339410,
        np.pi/180* 49.6198,
        np.pi/180* 1.8497,
        np.pi/180*(336.2075 - 49.6198),
        686.98,
        1.52371034*np.sqrt(1-np.power(0.09339410,2)),
        2*np.pi/686.98,
        2454942.90694,
        'red',
    ))

    #"""----------グラフ設定----------"""
    fig = plt.figure()
    ax = p3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    # fig.patch.set_facecolor('white')

    # 軸の設定
    ax.plot([-2,2], [0, 0], [0, 0] , color="g", lw=0.7)
    ax.plot([0,0], [-2, 2], [0, 0] , color="g", lw=0.7)
    ax.plot([0,0], [0, 0], [-2, 2] , color="g", lw=0.7)
    ax.set_xlim3d([-2,2])
    ax.set_ylim3d([-2,2])
    ax.set_zlim3d([-2,2])
    ax.view_init(elev=input_elev, azim=input_azim)

    # 惑星の軌道表示
    for planet in planets:
        t = np.linspace(0, 2*np.pi, 180)
        phase = 0
        xx = planet.a * np.cos(t - phase) - planet.a*planet.e
        yy = planet.b * np.sin(t - phase)
        x, y, z = coordtrans(xx,yy,planet.cosohm,planet.sinohm,planet.cosi,planet.sini,planet.cosomega,planet.sinomega)
        ax.plot(x, y, z, label='ellipse', c='black', lw=0.5)

    # 初期描画
    sun = ax.scatter(0, 0, 0, marker='.', c='red', s=150)
    drawing_planets = []
    for i, planet in enumerate(planets):
        x, y = ellipse(0, planet.a, planet.e)
        z = 0
        draw_obj = ax.plot(x,y,z, '.', c=planet.color)[0]
        drawing_planets.append(draw_obj)


    #"""----------初期値設定----------"""
    # 開始年月日
    Year1 = 2018
    Month1 = 7
    Day1 = 31
    Hour1 = 0
    Minute1 = 0
    JD1 = DAYtoJD(Year1,Month1,Day1,Hour1,Minute1)

    # 終了年月日
    Year2 = 2020
    Month2 = 10
    Day2 = 6
    Hour2 = 0
    Minute2 = 0
    JD2 = DAYtoJD(Year2,Month2,Day2,Hour2,Minute2)

    # 表示用テキスト
    YMD1 = 'START: ' + str(Year1) + '.' + str(Month1) + '.' + str(Day1)
    YMD2 = 'FINISH: ' + str(Year2) + '.' + str(Month2) + '.' + str(Day2)

    #"""----------メイン関数----------"""
    def update(step):
        JD = JD1 + step
        tt = JD
        for i, planet in enumerate(planets):
        # 地球
            u0 = planet.n*(tt-planet.t)
            u1 = u0
            for nn in range(10):
                u2 = u1-funct(u1,tt,planet.e,planet.n,planet.t)/diff_funct(u1,planet.e)
                u1 = u2
            xx, yy = ellipse_Kepller(u1, planet.a, planet.e)
            x, y, z = coordtrans(xx,yy,planet.cosohm,planet.sinohm,planet.cosi,planet.sini,planet.cosomega,planet.sinomega)
            drawing_planets[i].set_data(x, y)
            drawing_planets[i].set_3d_properties(z)
        return

    #"""---------グラフ設定(streamlit)----------"""
    # st.sidebar.warning('これより下作成中')
    # st.sidebar.text('惑星の表示')
    # Mercury =  st.sidebar.checkbox('水星')
    # Venus =  st.sidebar.checkbox('金星')
    # Earth = st.sidebar.checkbox('地球')
    # Mars = st.sidebar.checkbox('火星')
    # Jupiter = st.sidebar.checkbox('木星')

    #"""----------アニメショーンの実行----------"""
    # 実行フレームの設定
    frames = int(JD2 - JD1) + 1
    plot = st.pyplot(fig)
    st.write('{}'.format(YMD1),'　　', '{}'.format(YMD2))
    col1, col2, col3= st.columns(3)
    with col1:
        start = st.button('START')
    with col2:
        rerun = st.button('RERUN')
        if rerun:
            st.experimental_rerun()
    with col3:
        stop = st.button('STOP')
    
    if start:
        i = 0
        while i <= frames:
            update(i)
            plot.pyplot(fig)
            i += 3
        if stop:
            st.stop()

# r'''
# # ケプラー方程式による惑星会合周期の計算
# ## 火星の会合周期
# 火星と地球を太陽と中心とする円運動をすると仮定する。
# - 地球の公転周期:$E=365日$
# - 火星の公転周期:$P=687日$
# 角速度の差に会合周期Sをかけると、火星に対して地球の進んだ各360°となるから公転周期は以下のように計算できる。


# $(\frac{360°}{E} - \frac{360°}{P}) * S = 360°$

# $\frac{1}{E}  - \frac{1}{P} = \frac{1}{S}$

# $S = \frac{EP}{P - E} = \frac{365 * 687}{687 * 365} ≒ 779 ≒ 2年1.7ヶ月 ≒ 約2年2ヶ月$
# '''

main()
