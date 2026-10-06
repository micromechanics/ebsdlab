// piccolo_theme follows the operating system; show the light theme unless the reader has chosen one
try {
    if (!localStorage.getItem('piccoloThemeMode')) {
        localStorage.setItem('piccoloThemeMode', 'light');
    }
} catch (err) { /* storage blocked: theme falls back to its default */ }
