window.DEEPL_CONFIG = {
    targetLang: 'EN', // Default fallback target language (e.g., 'EN', 'ZH', 'JA', 'ES')
    // Elements/selectors that should NEVER be translated
    skipSelectors: [
        'code', 'pre', 'script', 'style', '.no-translate', 
        '.md-nav__link', '.md-footer', 'math', '.arithmatex'
    ],
    // DeepL API endpoint (Use 'https://api-free.deepl.com/v2/translate' if on the free plan)
    endpoint: 'https://api-free.deepl.com/v2/translate'
};