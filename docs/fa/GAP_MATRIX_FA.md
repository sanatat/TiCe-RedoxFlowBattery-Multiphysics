# پوشش گپ‌های پروژه Ti–Ce در نسخه جدید GitHub

این نسخه کدهای بیشتری برای نمایش توانایی **مدل‌سازی علمی** دارد؛ اما تجربه آزمایشگاهی جدید یا پیش‌بینی تأییدشده ایجاد نمی‌کند.

| گپ | توسعه جدید | وضعیت علمی |
|---|---|---|
| Mass Transfer | مدل انتقال دوبعدی قبلی + راکتور با مصرف Ce و Ti در امتداد جریان | حل عددی و تست موازنه؛ هندسه الکترود واقعی نیست |
| Reactor Engineering | افت فشار Darcy، توان پمپ، زمان اقامت، الکتروشیمی و نرخ مصرف واکنش‌دهنده | یک‌بعدی و با ضرایب سناریویی |
| CV/EIS | ماژول‌های قبلی با داده‌های جدول SI و منحنی‌های مصنوعی | هیچ آزمایش عملی انجام نشده |
| Ti–Ce Chemistry | سناریوهای ترکیب اسید MSA/H₂SO₄ و پیوند شرطی فلز–لیگاند | داده معتبر تعادل برای پیش‌بینی واقعی رسوب نداریم |
| Membrane & Crossover | نفوذ فیکی + حساسیت به میدان الکتریکی + جابه‌جایی آب، با حفظ موجودی عناصر | بار مؤثر گونه‌ها و شار آب اندازه‌گیری نشده |
| Battery Cycling | مدل قبلی + ۱۲۰۰ سیکل **فرضی** با تخریب برگشت‌پذیر و برگشت‌ناپذیر | ۳۷۲۱ سیکل واقعی مقاله را بازتولید نکرده‌ایم |
| CFD & Flow Fields | مدل دوبعدی قبلی + اثر دبی در راکتور و توزیع جریان استک | مدل هندسه واقعی یا CFD استک نیست |
| Scale-up & Stacking | مدل استک سری الکتریکی / موازی هیدرولیکی | جریان‌های نشتی الکتریکی و مانیفولد واقعی وارد نشده‌اند |

## اجرای توسعه‌ها

```bash
conda activate science
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python run_extensions.py
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/
```

خروجی‌ها در `extended_figures/` و `extended_outputs/` ساخته می‌شوند. نام و عنوان هر نمودار مشخص می‌کند «محاسبه شرطی»، «آزمایش مجازی» یا «داده SI» است.

**برای CV:** عبارت‌های علمی مناسب عبارت‌اند از *Finite-volume Mass Transfer, Coupled Electrochemical Reactor Model, Conditional Mixed-Acid Electrolyte Screening, Membrane Crossover Sensitivity, Stack Hydraulic Balancing, Long-Term Degradation Scenario Analysis*. تجربه «experimental CV/EIS» یا «validated precipitation prediction» ادعا نشود.
