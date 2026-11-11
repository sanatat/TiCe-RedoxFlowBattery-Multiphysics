# نسخه جدید انتقال جرم / Carbon Felt

این بسته شامل نسخه یکپارچه قبلی و لایه جدید برازش کاهش‌یافته شکل ۶c است. نمودارهای عددی جدید در `validation_figures/` و CSVها در `validation_outputs/` قرار می‌گیرند.

### Linux Mint

```bash
unzip ~/Downloads/TiCe-RedoxFlowBattery-Multiphysics.zip -d ~
cd ~/TiCe-RedoxFlowBattery-Multiphysics
conda activate science
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python calibrate_transport.py
OPENBLAS_NUM_THREADS=1 python design_tradeoffs.py
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/
```

اگر می‌خواهی مدل یکپارچه قبلی را هم از نو اجرا کنی:

```bash
OPENBLAS_NUM_THREADS=1 python run_integrated.py --fast
```

### چگونه داده‌های تقریباً استخراج‌شده را اصلاح کنم؟

در ابتدای فایل `calibrate_transport.py` دو آرایه `FIG6C_J` و `FIG6C_V` تعریف شده‌اند. این‌ها **اعداد ابزار اندازه‌گیری نیستند**، بلکه خواندن بصری از نمودار منتشرشده با خطای تخمینی ±۰٫۰۴ V هستند. بعد از هر اصلاح، اسکریپت را دوباره اجرا کن. نسخه اولیه داده‌ها در `validation_outputs/digitized_fig6c_approx.csv` هم ذخیره می‌شود.

### چرا سطح فعال Carbon Felt را اعلام نمی‌کنیم؟

برای تعیین آن لازم است حداقل هندسه واقعی میدان جریان، مقدار دقیق SOC هنگام اندازه‌گیری، مقاومت انتقال بار و یک اندازه‌گیری مستقل انتقال جرم یا سطح فعال را داشته باشیم. مدل فعلی **بازه سناریوهای سازگار** را نشان می‌دهد و ادعای اندازه‌گیری نمی‌کند.
