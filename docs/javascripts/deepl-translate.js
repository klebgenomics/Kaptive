async function translatePage(targetLang) {
    if (!window.DEEPL_API_KEY) {
        console.error("DeepL API Key is missing.");
        return;
    }

    const config = window.DEEPL_CONFIG;
    let textNodes = [];
    
    // 1. Walk the DOM and collect text nodes containing actual content
    const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT, {
        acceptNode: function(node) {
            // Check if node is inside a skipped selector
            for (let selector of config.skipSelectors) {
                if (node.parentElement && node.parentElement.closest(selector)) {
                    return NodeFilter.FILTER_REJECT;
                }
            }
            // Skip empty whitespace
            if (!node.nodeValue.trim()) return NodeFilter.FILTER_SKIP;
            return NodeFilter.FILTER_ACCEPT;
        }
    });

    while (walker.nextNode()) textNodes.push(walker.currentNode);
    if (textNodes.length === 0) return;

    // 2. Prepare texts for DeepL API payload
    const textValues = textNodes.map(node => node.nodeValue);

    // Build URLSearchParams as required by DeepL POST API
    const params = new URLSearchParams();
    params.append('target_lang', targetLang);
    textValues.forEach(text => params.append('text', text));

    try {
        // 3. Request translation from DeepL
        const response = await fetch(config.endpoint, {
            method: 'POST',
            headers: {
                'Authorization': `DeepL-Auth-Key ${window.DEEPL_API_KEY}`,
                'Content-Type': 'application/x-www-form-urlencoded'
            },
            body: params
        });

        if (!response.ok) throw new Error(`DeepL API Error: ${response.statusText}`);
        
        const data = await response.json();
        
        // 4. Update the DOM with translated strings
        data.translations.forEach((translation, index) => {
            textNodes[index].nodeValue = translation.text;
        });

        // Store preference locally so user settings persist across pages
        localStorage.setItem('deepl_lang_pref', targetLang);

    } catch (error) {
        console.error("Translation failed:", error);
    }
}

// 5. Intercept language switcher menu clicks
document.addEventListener("DOMContentLoaded", () => {
    window.addEventListener("hashchange", () => {
        const hash = window.location.hash;
        if (hash.startsWith("#deepl-translate-")) {
            const lang = hash.replace("#deepl-translate-", "").toUpperCase();
            translatePage(lang);
        }
    });

    // Auto-translate on page load if preference exists
    const savedLang = localStorage.getItem('deepl_lang_pref');
    if (savedLang) {
        // Add a small delay to ensure page fully loads
        setTimeout(() => translatePage(savedLang), 300);
    }
});