<script setup lang="ts">
import { ref } from 'vue';
import { useRoute, useRouter } from 'vue-router';
import { useUserStore } from 'src/stores/user';
import { useI18n } from 'vue-i18n';
import { useQuasar } from 'quasar';

const userStore = useUserStore();
const router = useRouter();
const route = useRoute();
const { t } = useI18n();
const $q = useQuasar();

const username = ref<string | null>(null);
const password = ref<string | null>(null);

const error = ref(false);

const login = async () => {
  error.value = false;
  try {
    if (username.value && password.value) {
      await userStore.sessionLogin({
        username: username.value,
        password: password.value,
      });
    } else {
      $q.notify({
        type: 'warning',
        message: t('message.enter_username_password'),
      });
    }
  } catch (err) {
    console.error(err);
  }

  if (await userStore.checkAuthentication()) {   // userStore.authenticated
    $q.notify({
      type: 'positive',
      message: t('message.successfully_logged_in'),
    });
    let next = '/'

    await router.push({ path: next });
  } else {
    error.value = true;
  }
};
</script>

<template>
  <section
    class="tw-flex tw-items-center tw-justify-center tw-h-screen tw-py-26 tw-bg-white tw-relative tw-overflow-hidden"
  >
    <img
      class="tw-absolute tw-bottom-0 md:tw-top-0 tw-right-0 lg:tw-h-full md:tw-h-full tw-w-1/2 md:tw-w-1/3"
      src="../assets/pattern-two-smashes-indigo-light-right.svg"
      alt=""
    />
    <img
      class="tw-absolute tw-top-0 tw-left-0 lg:tw-h-full md:tw-h-full w-1/2 md:tw-w-2/3"
      src="../assets/pattern-dots-2-indigo-light-left.svg"
      alt=""
    />
    <div class="tw-container tw-px-4 tw-mx-auto tw-relative">
      <div class="tw-max-w-lg tw-mx-auto">
        <div class="tw-text-center tw-mb-8">
          <a class="tw-inline-block tw-mx-auto tw-mb-6" href="#">
            <img src="../assets/nexus_logo.png" alt="" width="180" />
          </a>
          <h2
            class="tw-text-3xl md:tw-text-4xl tw-font-extrabold tw-mb-2 tw-text-indigo-500"
          >
            Sign in
          </h2>
        </div>
        <div>
          <div class="tw-mb-6">
            <label class="tw-block tw-mb-2 tw-font-extrabold" for="">{{
              t('label.username')
            }}</label>

            <input
              v-model="username"
              class="tw-inline-block tw-w-full tw-p-4 tw-leading-6 tw-text-lg tw-font-extrabold tw-placeholder-indigo-900 tw-bg-white tw-shadow tw-border-2 tw-border-indigo-900 tw-rounded"
              type="text"
            />
          </div>
          <div class="tw-mb-6">
            <label class="tw-block tw-mb-2 tw-font-extrabold" for="">{{
              t('label.password')
            }}</label>
            <input
              v-model="password"
              class="tw-inline-block tw-w-full tw-p-4 tw-leading-6 tw-text-lg tw-font-extrabold tw-placeholder-indigo-900 tw-bg-white tw-shadow tw-border-2 tw-border-indigo-900 tw-rounded"
              type="password"
            />
          </div>

          <button
            @click="login"
            class="tw-inline-block tw-w-full tw-py-4 tw-px-6 tw-mb-6 tw-text-center tw-text-lg tw-leading-6 tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 shadow tw-rounded tw-transition tw-duration-200"
          >
            {{ t('label.login') }}
          </button>
        </div>
      </div>
    </div>
  </section>
</template>

<style scoped>
@import 'tailwindcss/base.css';
@import 'tailwindcss/components.css';
@import 'tailwindcss/utilities.css';
</style>
