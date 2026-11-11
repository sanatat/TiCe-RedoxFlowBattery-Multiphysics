# راهنمای انتشار در GitHub

این پوشه **یک Repository آماده** است. نیازی نیست PDF مقاله یا فایل Supporting Information را در GitHub آپلود کنی. خروجی‌ها و نمودارهای تولیدشده داخل پوشه هستند؛ تصویر کپی‌شده از شکل مقاله عمداً حذف شده است.

## مراحل ساده (روش وب)

۱. در GitHub یک مخزن با نام `TiCe-RedoxFlowBattery-Multiphysics` بساز. برای اینکه فایل‌ها مستقیم داخل ریشه قرار بگیرند، هنگام ساخت مخزن README یا License جدید اضافه نکن؛ این پروژه از قبل هر دو فایل را دارد.

۲. فایل ZIP را روی لپ‌تاپ استخراج کن و پوشه را با دستورهای زیر در GitHub منتشر کن. به‌جای `YOUR_USERNAME` نام کاربری خودت را بگذار.

```bash
cd ~/TiCe-RedoxFlowBattery-Multiphysics
git init
git add .
git commit -m "Initial release: Ti-Ce RFB multiphysics research prototype"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/TiCe-RedoxFlowBattery-Multiphysics.git
git push -u origin main
```

اگر Git روی سیستم هویت Commit را تنظیم نکرده باشد، از نام و ایمیلی که می‌خواهی به Commit متصل شود با `git config --global user.name` و `git config --global user.email` استفاده کن. برای HTTPS ممکن است ورود از طریق مرورگر یا یک روش احراز هویت سازگار با GitHub لازم باشد؛ از پسورد حساب به جای توکن استفاده نکن.

## عنوان و توضیح پیشنهادی

**Repository name:** `TiCe-RedoxFlowBattery-Multiphysics`

**Description:** متن کوتاه آماده‌شده در فایل [`GITHUB_METADATA.md`](GITHUB_METADATA.md) را در بخش About قرار بده؛ این متن برای محدودیت طول توضیح GitHub تنظیم شده است.

**Topics:** فهرست پیشنهادی در [`GITHUB_METADATA.md`](GITHUB_METADATA.md) موجود است؛ روی موضوعات واقعاً شبیه‌سازی‌شده تمرکز دارد.

## تست قبل از Push

```bash
python -m pip install -r requirements.txt
python run_all.py
python -m pytest -q tests/
```

به‌صراحت در رزومه عنوان **Computational Research Project** استفاده کن، نه «Experimental RFB Research».

## توسعه‌های جدید

نسخه به‌روزشده شامل `extensions/`، `run_extensions.py`، `extended_figures/`، `extended_outputs/`، `TiCe_Gap_Extensions.ipynb` و `MODEL_COVERAGE.md` است. برای بررسی مستقل:

```bash
python run_extensions.py
python -m pytest -q tests/
```

**۶۶ تست عددی** باید پاس شوند. آمار تست‌ها و خروجی‌ها مربوط به محیط اجرای این نسخه است. چگالی انرژی واقعی، شیمی رسوب Ti/Ce، بازده سلول آزمایشگاهی و طول عمر واقعی به داده‌های مستقل نیاز دارند؛ نتایج سناریویی را در GitHub به‌عنوان آزمایش گزارش نکن.

## افزونه‌ی جدید غربالگری ترمودینامیکی

از منابع منتشرشده برای حلالیت Ce(III)/Ce(IV) در ترکیب‌های مشخص MSA/H2SO4 و یک مدل تعادلی محدود به محیط 2 M NaClO4 استفاده شده است. **مدل رسوب قطعی Ti/Ce در اسید مخلوط نیست.** دستور اجرای جداگانه:

```bash
python run_evidence_thermodynamics.py
python -m pytest -q tests/
```

فایل `THERMODYNAMICS_EVIDENCE.md` مراجع و محدودیت‌های علمی را فهرست کرده است.


## README جدید و کلیدواژه‌های آگهی

متن انگلیسی `README.md` به‌صورت طبیعی و پژوهشی بازنویسی شده است. کلیدواژه‌های مهم مانند **Mass Transfer، Reactor Engineering، Titanium–Cerium Redox Chemistry، Electrolyte Formulation، Flow Field Design، Membrane Crossover، Power Density، Battery Cycling، Scale-up، Stacking و Long-Term Stability** در توضیح کاربرد واقعی مدل‌ها آمده‌اند، نه صرفاً یک لیست تکراری برای جست‌وجو.

فایل [`RESEARCH_KEYWORD_MAP.md`](RESEARCH_KEYWORD_MAP.md) ارتباط هر کلیدواژه با فایل کد و خروجی را روشن می‌کند. بخش‌هایی که واقعاً انجام نشده‌اند، از جمله **آزمایش عملی CV/EIS، ساخت فیزیکی الکترود، و پیش‌بینی معتبر رسوب Ti/Ce** به‌صورت صریح از توانایی‌های پیاده‌سازی‌شده جدا شده‌اند.

## نسخه انتشار

بعد از استخراج ZIP نهایی، اجرای `python run_all.py` و `python -m pytest -q tests/` را بررسی کن و سپس پوشه را Push کن. لازم نیست PDF مقاله یا SI را داخل مخزن آپلود کنی.
