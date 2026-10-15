function onEnter(e, submit) {
    if (e.which == 13) {
        console.log("clicked " + submit);
        document.getElementById(submit).click();
    }
}

function handleEvent(submit) {
    return function(e) {
        onEnter(e, submit);
    };
}

function bindEnter(id, submit) {
    document.getElementById(id).addEventListener('keypress', handleEvent(submit));
}