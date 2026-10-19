/**
 * Disable all buttons while Shiny is busy.
 *
 * Listens for Shiny lifecycle events:
 *  - 'shiny:busy': triggered when the server is processing
 *  - 'shiny:idle': triggered when processing completes
 *
 * Behaviour:
 *  - Disables all <button> elements while busy
 *  - Re-enables them when idle
 *
 * Purpose:
 *  - Prevent duplicate clicks / repeated submissions
 *  - Avoid overlapping reactive execution
 *
 * Caveats:
 *  - Applies globally to all buttons
 *  - May disable unrelated UI controls
 */

$(document).on("shiny:busy", function() {
  var inputs = document.getElementsByTagName("button");
  console.log(inputs);
for (var i = 0; i < inputs.length; i++) {
inputs[i].disabled = true;
}
});

$(document).on("shiny:idle", function() {
  var inputs = document.getElementsByTagName("button");
  console.log(inputs);
for (var i = 0; i < inputs.length; i++) {
inputs[i].disabled = false;
}
});
