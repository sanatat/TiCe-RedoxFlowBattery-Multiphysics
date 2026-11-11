# NEW: Transport Calibration + Identifiability Audit

**Start here:** `QUICKSTART_TRANSPORT_FA.md` and `REPORT_MT_CALIBRATION_FA.md`.

New package adds discharge-polarization Figure 6c inverse analysis (approximate manual digitization), 18 finite-volume Ti/Ce flow-diffusion scenarios, conditional flow/pump design comparisons, structural porous-felt sensitivity, held-out internal assessment and 12 extra tests (total **37**). No physical active-surface measurement or independently validated mass-transfer limit is claimed.

---

# Ti–Ce Redox Flow Battery — Integrated Transport & Electrochemistry Study

**بازسازی محاسباتی آموزشی از مقاله Fu et al. (2026) و فایل Supplementary (۴۱ صفحه)**

**DOI مقاله:** https://doi.org/10.1002/cssc.202502732

**وضعیت علمی:** پروژه عددی قابل اجرا و دارای تست‌های محافظت موازنه جرم است؛ **Digital Twin اعتبارسنجی‌شده یا کار آزمایشگاهی نیست.** محاسبات جدید CFD، انتقال جرم، CV مجازی، EIS مجازی و سیکل باتری روی سخت‌افزار محلی انجام می‌شوند. کالیبراسیون واقعی تمام پارامترها و بازتولید ۳۷۲۱ سیکل مقاله انجام نشده است.

## اجرای پروژه (Linux Mint)

```bash
unzip TiCe-RedoxFlowBattery-Multiphysics.zip
cd TiCe_RFB_Integrated
conda activate science
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python run_integrated.py
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/
```

برای اجرای سریع روی CPU قدیمی:

```bash
OPENBLAS_NUM_THREADS=1 python run_integrated.py --fast
```

نتایج در `figures/` و `outputs/` ذخیره می‌شوند. نتیجه قابل مشاهده در `outputs/summary.json` است.

## پوشش هر هفت مهارت

| حوزه | مدلی که واقعاً اجرا می‌شود | فایل کد | محدودیت باقی‌مانده |
|---|---|---|---|
| Mass Transfer | حل دوبعدی **Advection–Diffusion–Reaction** با روش Finite Volume، شار دیواره، حد جریان و موازنه دقیق جرم؛ تحلیل اثر D و Q | `transport_fv.py`, `run_integrated.py` | هندسه الکترود واقعی معلوم نیست؛ D گزارش‌شده ظاهری است |
| Reactor Engineering | Re، فشار، توان پمپ، زمان اقامت، گردش مخزن و تبدیل گونه‌ها | `transport_fv.py` | کانال و راندمان پمپ فرضی‌اند؛ افت فشار واقعی Carbon Felt اندازه‌گیری نشده |
| Electrochemistry | CV سه‌الکترودی مجازی با نفوذ و BV، رگرسیون جدول تجربی S4/S5، EIS مدار معادل مجازی | `cv_virtual.py`, `electrochem_cell.py` | **CV/EIS واقعی انجام نشده**؛ پارامتر Cdl و Rcontact فرضی است |
| Ti–Ce Electrolyte Chemistry | موازنه‌های Ce(III)/Ce(IV) و Ti(III)/Ti(IV)، Nernst و داده‌های واقعی رسانایی/ویسکوزیته در غلظت و دما | `electrochem_cell.py` | **مدل گونه‌بندی، پایداری رسوب و اسیدهای مخلوط اعتبارسنجی نشده**؛ DFT تکرار نشده |
| Membrane & Crossover | مدل Fick با نفوذپذیری گزارش‌شده Ce، موازنه موجودی Ce/Ti، سناریوی Ti مجهول، بازتنظیم ایده‌آل | `electrochem_cell.py` | هندسه H-cell، انتقال پروتون و آب، و نفوذپذیری دقیق Ti نامشخص‌اند |
| Battery Cycling | واکنش فارادایی، Butler–Volmer، Nernst، افت اهمی، افت غلظتی محدود به نتایج انتقال جرم دوبعدی، بازده CE/VE/EE، سناریوی بازیابی | `electrochem_cell.py` | ۱۲ سیکل سناریویی، نه آزمایش ۳۷۲۱ سیکلی؛ سینتیک تخریب واقعی مشخص نیست |
| CFD & Flow Fields | D2Q9 LBM سرعت دوبعدی با/بدون مانع + حل FVM انتقال یون با همان میدان سرعت | `lbm_flow.py`, `transport_fv.py` | هندسه فرضی پریودیک است؛ نسبت مقاومت هیدرولیکی **مقیاس‌گذاری تقریبی** است |

### اضافه: مدل واکنش/نفوذ داخل Carbon Felt

در `porous_felt.py` از ضخامت **واقعی** ۵ میلی‌متری SI استفاده شده است؛ اما **تخلخل، tortuosity و سطح فعال داخلی فرضی** هستند. عدد تیله، ضریب اثربخشی و عمق نفوذ واکنش حساب می‌شود. حل تحلیلی و حل Finite Volume مستقل با هم مقایسه شده‌اند. سپس ضریب موثر واکنش در مرز، به حل دوبعدی کانال متصل می‌شود.

## واقعی یا فرضی؟

**مستقیم از SI:** مساحت ۹cm²، حجم هر مخزن ۲۰mL، دبی هر سمت ۲۰mL/min، Carbon Felt با ضخامت ۵mm، غشای N115، ۲M MSA، ۰.۵M هر گونه فعال، ولتاژ قطع ۱.۸ و ۰.۸V، D و k ظاهری جدول S6، رسانایی/ویسکوزیته S7/S8، نفوذپذیری Ce از شکل S4 و مقاومت غشا از شکل S20.

**فرضی/سناریویی:** ابعاد کانال ۳۰×۱×۳۰mm، سطح ویژه الکترود متخلخل، تخلخل و tortuosity، بازده پمپ، فاصله عبور یونی، مقاومت تماسی، ظرفیت دو‌لایه، نرخ واکنش شیمیایی هنگام اختلاط، ضریب عبور Ti و هندسه دقیق جریان. هیچ‌کدام را به مقاله نسبت نمی‌دهیم.

**ناسازگاری منبع:** متن SI صفحه ۳ ترکیب mRFB را مانند sRFB تکرار کرده؛ با متن مقاله اصلی و شکل S10 مغایرت دارد. در کد، ترکیب mRFB طبق توضیح مقاله اصلی و SI S10 در هر دو مخزن **۰.۵M Ce + ۰.۵M Ti + ۲M MSA** فرض شده و این تصمیم مستند است.

## نکات مهم علمی

1. ضریب‌های D از CV الکترود متخلخل و واکنش شبه‌برگشت‌پذیر به دست آمده‌اند. **نباید** آن‌ها را بدون اعتبارسنجی ضرایب نفوذ حجمی واقعی فرض کرد. نمودار حساسیت D×۱، D×۰٫۱ و D×۰٫۰۱ آماده شده است.
2. در شرایط انتقال جرم، فرض دیواره کامل جذب‌کننده **کران بالای انتقال** است. مدل Thiele سناریوی محدودیت واکنش داخل Felt را نشان می‌دهد.
3. EIS کاملاً **مجازی و مبتنی بر مدار معادل** است؛ این فایل داده آزمایشگاهی EIS ندارد.
4. فقط دو مقدار دقیق EE در متن مقاله برای rate-test داریم: 86.02% در 30 و 64.37% در 90 mA/cm². در شکل ۱۳ مقایسه شده‌اند؛ تطابق کامل مشاهده نشده است و **مدل با آن‌ها fit نشده**.
5. در شبیه‌سازی Regeneration، نفوذپذیری Ce عمداً ۳۰۰ برابر زیاد شده است تا اثر بازیابی در ۱۲ سیکل قابل مشاهده شود. این یک **Stress Test فرضی** است، نه ادعای نتیجه واقعی مقاله.
6. مدل پس از اختلاط، واکنش Ce(IV) + Ti(III) را به‌صورت استوکیومتری و فوری اجرا می‌کند؛ در مقاله ثابت سرعت این واکنش موجود نیست.

برای توضیح حدود اعتبارسنجی به `VALIDATION_AND_GAPS_FA.md` و برای کل فهرست شکل‌ها به `FIGURE_INDEX.md` مراجعه کن.

### داده‌های محاسبات کوانتومی مقاله

در `electrolyte_chem.py` مقادیر **الکترونیکی** انرژی پیوند کمپلکس‌های Ce(III)، Ce(IV)، Ti(III) و ساختار دوهسته‌ای Ti(IV) از جدول S2 وارد شده است. نمودار `14_reported_DFT_coordination_energies.png` این داده‌ها را ترسیم می‌کند. این‌ها انرژی آزاد گیبس یا ثابت تعادل نیستند و **نمی‌توان** از آن‌ها به‌تنهایی غلظت گونه‌ها یا رسوب را پیش‌بینی کرد. علاوه‌بر این، رابطه Nernst برای جفت Ti(IV)/Ti(III) شامل دو پروتون به‌صورت **حساسیت ایده‌آل نسبت به فعالیت H+** نوشته شده است؛ غلظت اسید MSA الزاماً برابر فعالیت H+ نیست.
