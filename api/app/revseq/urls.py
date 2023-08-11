from django.contrib import admin
from django.urls import path, include, re_path
from users.views import CsrfCookieView, LoginView, LogoutView, UserViewSet
from core.views import (
    SampleCountViewSet,
    PlateListView,
    MetadataViewSet,
    SubstrainListView,
)
from rest_framework.routers import DefaultRouter

router = DefaultRouter()
router.register(r"api/samplecounts", SampleCountViewSet, basename="samplecounts")
router.register(r"api/metadata", MetadataViewSet, basename="metadata")

urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/auth/cookie/", CsrfCookieView.as_view(), name="auth-cookie"),
    path("api/auth/login/", LoginView.as_view(), name="login"),
    path("api/auth/logout/", LogoutView.as_view(), name="logout"),
    path("api/auth/users/me/", UserViewSet.as_view({"get": "retrieve"}), name="me"),
    path("api/plates/", PlateListView.as_view(), name="plates"),
    path("api/substrains/", SubstrainListView.as_view(), name="substrains"),
    path("", include(router.urls)),
]
