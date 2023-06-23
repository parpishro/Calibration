
#
#
# X <- data.frame(height=ballField[,2], time=ballField[,1])
# X2 <- data.frame(height=ballSim[,2], time=ballSim[,1])
# X2 <- X2[order(X2$height), ]
#
# plot(X[, 1], X[, 2], cex = 0.65, pch = 19, col = "red", ylim = c(0, 1.5),
#      xlab="Height (m)", ylab="Time (s)")
# # lines(X[1:21, 1], X[1:21, 3], cex = 0.65, pch = 19, col = "blue")
# # lines(X[1:21, 1], X[1:21, 3]+(2*X[1:21, 4]), cex = 0.65, pch = 19, col = "grey")
# # lines(X[1:21, 1], X[1:21, 3]-(2*X[1:21, 4]), cex = 0.65, pch = 19, col = "grey")
# lines(X[1:21, 1], preds$pred[1:21], cex = 0.65, pch = 19, col = "blue")
# lines(X[1:21, 1], preds$pred[1:21]+(2*preds$se[1:21]), cex = 0.65, pch = 19, col = "grey")
# lines(X[1:21, 1], preds$pred[1:21]-(2*preds$se[1:21]), cex = 0.65, pch = 19, col = "grey")
# lines(X2[, 1], X2[, 2], cex = 0.65, pch = 19, col = "green")
#
#
# legend("topleft", legend = c("Field", "Prediction", "PI", "Simulation"),
#        col = c("red", "blue", "grey",  "green"), pch = 16, cex = 0.5)
