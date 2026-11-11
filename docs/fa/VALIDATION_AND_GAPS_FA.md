# گزارش اعتبار عددی و گپ‌های باقی‌مانده

## موفقیت‌های محاسباتی

- روش Finite Volume طوری نوشته شده که شار مولی ورودی = شار خروجی + شار واکنش دیواره باشد. اختلاف عددی در مرتبه شناوری دستگاه است.
- مدل انتقال از تغییر دبی، ضریب نفوذ و واکنش دیواره پاسخ دارد و با مش ظریف‌تر از آزمون همگرایی عبور می‌کند.
- مدل ۱D درون Felt با راه‌حل تحلیلی واکنش–نفوذ (Thiele/effectiveness factor) و حل مستقل ماتریسی مقایسه می‌شود.
- مدل LBM سرعت را برای کانال ساده و دارای مانع حل می‌کند؛ حلال انتقال گونه از میدان سرعت تولیدی استفاده می‌کند. موازنه جرم جزء منتقل‌شده بررسی می‌شود.
- واکنش Faradaic و انتقال غشایی هرکدام موجودی کل Ti و Ce را حفظ می‌کنند. اختلاط ایده‌آل هم همین شرط را حفظ می‌کند.
- نمودارهای CV مجازی و EIS مجازی تولید می‌شوند؛ جدول پیک‌های CV *واقعی مقاله* S4/S5 جداگانه تحلیل می‌شود.
- اثر غلظت و دما بر رسانایی/ویسکوزیته از جدول واقعی S7/S8 آورده می‌شود.

## چیزهایی که نمی‌شود بدون داده جدید ادعا کرد

| گپ | دلیل | داده لازم برای تکمیل |
|---|---|---|
| تطبیق Mass Transfer با الکترود واقعی | سطح داخلی/تخلخل/تورتوزیته و ساختار flow-field دقیق در SI نیست | تصویربرداری میکروساختار، BET مناسب، porosimetry یا برآورد سطح موثر از آزمایش limiting current |
| کالیبراسیون سینتیک الکترود | k از CV در Carbon Felt به سطح هندسی وابسته است؛ Rct واقعی معلوم نیست | CV خام/EIS خام با هندسه شناخته‌شده و iR compensation |
| CFD واقعی استک | هندسه کانال، توزیع سیال و افت فشار معلوم نیست | CAD/مش واقعی، فشار تفاضلی و منحنی Δp–Q |
| گونه‌بندی Ti/Ce در MSA و Mixed-acid | ثابت‌های تعادل و activity model معتبر موجود نیست | ثابت‌های تشکیل کمپلکس، داده‌های حلالیت/رسوب در اسید و دما |
| طول عمر و Regeneration | Crossover تابع SOC و میدان الکتریکی/آب؛ واکنش‌های جانبی و تخریب Felt مشخص نیست | داده‌های ۳۷۲۱ سیکل، حجم/غلظت هر مخزن بر حسب زمان و تاریخچه Mixing |
| EIS واقعی | خود مقاله منحنی EIS کامل ارائه نکرده؛ فقط ASR غشا | طیف EIS سلول واقعی در چند SOC و دبی |
| مقایسه معتبر EE | EE تابع پروتکل، پیش‌سیکل، SOC اولیه و شرایط ناشناخته است | ولتاژ-جریان-زمان خام، روش دقیق پیش‌شرطی‌سازی و نرخ‌های شارژ/دشارژ |

## فرمول‌های مدل

**انتقال جرم:**

`div(u*c - D*grad(c)) = 0`, در ورودی `c=cin` و در دیواره واکنش‌پذیر `-n·(u*c-D*grad(c)) = k_surface*c_surface`.

**مدل Felt:** `D_eff = D*epsilon/tortuosity`, `phi = L*sqrt(a_s*k_s / D_eff)`, `effectiveness = tanh(phi)/phi`, `k_wall = sqrt(D_eff*a_s*k_s)*tanh(phi)`.

**انتقال غشا:** `dn_i/dt = ± P_i*A_m*(c_positive-c_negative)` برای هر گونه. این مدل به‌تنهایی انتقال الکترومهاجرتی پروتون و تغییر حجم را ندارد.

**راکتور:** `Re = rho*u*(2H)/mu` و `delta_p = 12*mu*L*u/H²` برای کانال صفحه‌ای ایده‌آل (نه Carbon Felt واقعی).

**ولتاژ:** اختلاف پتانسیل‌های دو نیم‌سلول از Nernst و Butler–Volmer معکوس، به‌علاوه افت‌های Ohmic و Concentration overpotential.

**داده SI:** diffusion و kinetics جدول S6، Conductivity/viscosity جدول S7/S8، N115 و H-cell شکل S4/S20. همه تبدیل‌های واحد SI در کد انجام شده‌اند.

## دقت علمی و کاربری برای CV

عنوان قابل دفاع پس از اجرای شخصی کد:

**Independent computational study — Ti–Ce redox flow battery transport and performance (Python)**

عبارت قابل دفاع: *Implemented and tested conservative 2D advection–diffusion–reaction models, an idealized porous-felt reaction–diffusion submodel, membrane crossover balances, virtual electrochemical diagnostics, and a coupled charge–discharge performance simulator using published supplemental parameters; documented model identifiability and validation limitations.*

**ادعا نکن:** experimental CV/EIS, battery fabrication, 3,721-cycle reproduction, verified pore-scale CFD, validated mixed-acid Ti–Ce speciation, or a true digital twin.
