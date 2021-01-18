
function cpp_output() {
    /** @type {HTMLCanvasElement} */
    var canvas = document.getElementById("canvas-cpp-output");
    var ctx = canvas.getContext("2d");

    update_canvas_margins(canvas);
    window.addEventListener("resize", function() { update_canvas_margins(canvas); });

    /** @type {CanvasRenderingContext2D} */
    ctx.fillStyle = "rgb(255, 255, 255)";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
}

function update_canvas_margins(canvas) {
    canvas.width = 31*window.innerWidth/32;
    canvas.height = 31*window.innerHeight/32;
}
