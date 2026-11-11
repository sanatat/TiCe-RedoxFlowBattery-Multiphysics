from pathlib import Path
import nbformat as nbf
root=Path(__file__).resolve().parent
nb=nbf.v4.new_notebook()
nb.metadata={'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'},'language_info':{'name':'python'}}
cells=[]
cells.append(nbf.v4.new_markdown_cell('''# Ti–Ce RFB — Integrated Battery Modeling\n\nمدل محاسباتی آموزشی بر اساس مقاله Fu et al. 2026 و SI؛ **نه آزمایشگاه و نه Digital Twin اعتبارسنجی‌شده**.\n\n**پوشش:** CFD و انتقال جرم دوبعدی، واکنش داخل Felt، مهندسی راکتور، CV/EIS مجازی، شیمی الکترولیت، Crossover، سیکل باتری.\n\nمرجع: https://doi.org/10.1002/cssc.202502732\n\nبرای اجرای همه آزمایش‌ها با نمودارهای قابل تکرار، `python run_integrated.py` را اجرا کن. این Notebook مسیر اصلی را مرحله‌به‌مرحله بررسی می‌کند.'''))
cells.append(nbf.v4.new_code_cell('''import json, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from transport_fv import Channel, channel_transport, lbm_velocity
from porous_felt import Felt, effectiveness, finite_difference_flux
from electrochem_cell import (Battery, SI, mass_transport_inputs, run_cycles,
    virtual_eis, electrolyte_props, published_apparent_D)
print('DOI:', SI['citation']['doi'])
print('SI geometry:', SI['cell']['area_cm2'], 'cm²,', SI['cell']['reservoir_each_mL'], 'mL each,',
      SI['cell']['flow_each_mL_min'], 'mL/min each')'''))
cells.append(nbf.v4.new_markdown_cell('''## ۱. انتقال جرم (محور اصلی)\n\nحل محافظه‌کار دوبعدی با معادله `div(uc − D grad(c)) = 0` و شرط مرزی واکنش روی الکترود. D مربوط به جدول S6 و **ظاهری** است.'''))
cells.append(nbf.v4.new_code_cell('''ch=Channel()
D=published_apparent_D('Ti3')
t=channel_transport(ch,D,250.,nx=90,ny=32)
print('Inlet, Outlet, reacted mol/s:', t['in_mol_s'],t['out_mol_s'],t['reacted_mol_s'])
print('Relative molar balance error:', t['mass_balance_relative_error'])
plt.figure(figsize=(8,3.5))
plt.pcolormesh(t['x_m']*1000,t['y_m']*1000,t['c'],shading='nearest')
plt.colorbar(label='Ti(III) mol/m³');plt.xlabel('x (mm)');plt.ylabel('y (mm)');plt.title('2-D conservative transport');plt.show()'''))
cells.append(nbf.v4.new_markdown_cell('''## ۲. انتقال جرم در الکترود متخلخل ۵ میلی‌متری\n\nتخلخل، tortuosity و سطح مؤثر داخلی در SI گزارش نشده‌اند؛ مقادیر زیر صرفاً فرضی‌اند.'''))
cells.append(nbf.v4.new_code_cell('''for a in [2500.,25000.,100000.]:
    f=Felt(active_area_m2_m3=a)
    d=effectiveness('Ti3',f)
    independent=finite_difference_flux('Ti3',f,n=500)
    linked=channel_transport(ch,D,250.,k_wall_m_s=d['k_wall_m_s'],nx=90,ny=32)
    print(f'a_s={a:.0f} m²/m³, phi={d["Thiele"]:.2f}, effectiveness={d["effectiveness"]:.4f},'
          f' FVM relative error={independent["relative_error"]:.3%},'
          f' cell equivalent j={linked["j_equiv_mA_cm2"]:.1f} mA/cm²')'''))
cells.append(nbf.v4.new_markdown_cell('''## ۳. مهندسی راکتور و افت فشار\n\nابعاد کانال فرضی‌اند، حجم و دبی مخازن از SI هستند. افت فشار مربوط به مدل کانال صفحه‌ای است، نه به کل Felt.'''))
cells.append(nbf.v4.new_code_cell('''print('Reynolds:', ch.reynolds)
print('Channel pressure loss (Pa):',ch.pressure_drop_Pa)
print('Pump power two cells (W):',ch.pump_power_two_cells_W)
print('Tank turnover (s):',Battery().tank_m3/Battery().flow_m3_s)
for flow in [10,20,40]:
    pp=Channel(flow_m3_s=flow*1e-6/60)
    rr=channel_transport(pp,D,250,nx=75,ny=28)
    print('Q',flow,'mL/min','conversion',rr['reacted_mol_s']/rr['in_mol_s'], 'k',rr['k_apparent_m_s'])'''))
cells.append(nbf.v4.new_markdown_cell('''## ۴. اتصال انتقال جرم به باتری\n\nحد جریان در هر نیم‌سلول از حل دوبعدی می‌آید؛ سینتیک ولتاژ از مدل Butler–Volmer و خواص SI است. ریزساختار Carbon Felt واقعی هنوز نامعلوم است.'''))
cells.append(nbf.v4.new_code_cell('''b=Battery()
km, _=mass_transport_inputs(b,nx=80,ny=30)
print('Coupled wall transfer coeffs (m/s):',km)
trace,summ,final,audit=run_cycles(b,km,cycles=3,j_mA_cm2=50,regenerate_after=(2,))
for row in summ: print('cycle',row['cycle'],'CE',round(row['CE_percent'],2),'EE',round(row['EE_percent'],2),
                       'capacity mAh',round(row['discharge_Ah']*1000,2))
print('Ce/Ti inventory conserved:', audit['initial_inventory_mol'],audit['final_inventory_mol'])'''))
cells.append(nbf.v4.new_markdown_cell('''## ۵. CV و EIS مجازی\n\nجدول‌های S4/S5 داده CV واقعی مقاله هستند؛ EIS زیر فقط مدار معادل محاسباتی است و نباید به آزمایشگاه نسبت داده شود.'''))
cells.append(nbf.v4.new_code_cell('''z=virtual_eis(b,km)
plt.figure(figsize=(6,3.5))
zz=z['Z_area_ohm_m2']*1e4
plt.plot(zz.real,-zz.imag);plt.xlabel("Z' (Ω cm²)");plt.ylabel("-Z'' (Ω cm²)");
plt.title('Virtual EIS only');plt.show()
print('SI CV table S4 first row:', {k:v[0] for k,v in SI['CV_Ce'].items() if isinstance(v,list)})'''))
cells.append(nbf.v4.new_markdown_cell('''## ۶. شیمی الکترولیت و محدودیت اساسی\n\nرسانایی و ویسکوزیته از S7/S8 استخراج شده‌اند. **گونه‌بندی مولکولی واقعی، حلالیت رسوب و mixed-acid thermodynamics مدل نشده‌اند**؛ باید ثابت‌های تشکیل کمپلکس و داده‌های تجربی مستقل داشته باشیم.'''))
cells.append(nbf.v4.new_code_cell('''for M in [.5,.75,1.]:
    print('Each metal',M,'M, sigma(S/m),mu(Pa s):',electrolyte_props(M))'''))
cells.append(nbf.v4.new_markdown_cell('''## ۷. CFD + FVM مشترک\n\nحل LBM کانال فرضی با مانع و انتقال Ti(III) تحت میدان سرعت حل‌شده؛ دوره‌ای و با نگاشت هندسه فرضی است، نه CAD واقعی استک.'''))
cells.append(nbf.v4.new_code_cell('''velx,vely,solid,flow=lbm_velocity(ch,nx=75,ny=27,baffles=True,steps=1200)
tr=channel_transport(ch,D,250.,flow_field=(velx,vely,solid))
print('Conservation error under LBM velocity:',tr['mass_balance_relative_error'])
plt.figure(figsize=(7,3.2));plt.imshow(np.ma.masked_array(tr['c'],mask=solid),origin='lower',aspect='auto')
plt.colorbar(label='Ti(III) mol/m³');plt.title('LBM/FV transport through illustrative baffles');plt.show()'''))
cells.append(nbf.v4.new_markdown_cell('''## نتیجه\n\nاین Notebook شواهد **توانایی مدل‌سازی محاسباتی** را نشان می‌دهد؛ نه تجربه عملی CV، EIS، مونتاژ باتری یا بازسازی دقیق ۳۷۲۱ سیکل. برای گزارش علمی و گپ‌های نیازمند داده جدید، `README_FA.md` و `VALIDATION_AND_GAPS_FA.md` را بخوان.'''))
nb['cells']=cells
nbf.write(nb,root/'TiCe_Integrated_Study.ipynb')
print(root/'TiCe_Integrated_Study.ipynb')
